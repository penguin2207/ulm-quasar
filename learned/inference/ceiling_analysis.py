#!/usr/bin/env python
# SPDX-License-Identifier: MIT
"""Ceiling-breakpoint + linear-range-end + bootstrap-CI analysis for the cross-algorithm count.

Per dataset (jun23 = ceiling test, apr17 = negative control), method (classical/deepulm/lista) and
tube (jun23: combinedTube/tubeL/tubeR; apr17: combinedTube):
  - per-window log-log slope + 95% CI
  - Muggeo-style continuous two-segment breakpoint (grid-searched) + bootstrap 95% CI; slope below/above
  - LINEAR-RANGE-END: highest concentration whose local (sliding-window) slope is still consistent with
    >= linear (upper CI bound >= 1, i.e. has not rolled below proportional). Per-resample criterion =
    sliding-window slope >= 1; reported with bootstrap 95% CI. This is the "where does the count leave
    the linear/quantitative range" ceiling, which is cleaner than the above-plateau slope and is robust to
    the learned methods being SUPRALINEAR (slope > 1) at low-mid concentration.

Conventions match the classical RUNNER exactly (verified): in-tube localization RATE, Bg-pedestal-
subtracted (mean over bgFlow blocks), no tracking/QC, BLOCK-LEVEL bootstrap (resample blocks within each
rung AND the bg blocks). jun23 excludes M3 (rung idx 7). Classical per-block from readout.locRate/bgLoc;
learned from count_locrate per_block (fixed pre-registered Bg-calibrated SR-map threshold).
"""
import numpy as np, h5py, os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import datasets as DS
from paths import OUTPUT_ROOT, SCRATCH_DIR  # machine-local roots; see learned/README.md

NBOOT = 1500
SEED = 1

# per-dataset analysis config
DSC = {
    "jun23": dict(exclude=[7], slide_w=5, tubes=("combinedTube", "tubeL", "tubeR"),
                  windows={"L1-L5": [0, 1, 2, 3, 4], "M1-M5": [5, 6, 8, 9],
                           "U1-U5": [10, 11, 12, 13, 14], "all": [i for i in range(15) if i != 7]}),
    "apr17": dict(exclude=[], slide_w=5, tubes=("combinedTube",),
                  windows={"C1-C4": [0, 1, 2, 3], "C5-C8": [4, 5, 6, 7],
                           "C3-C8": [2, 3, 4, 5, 6, 7], "all": [0, 1, 2, 3, 4, 5, 6, 7]}),
}


def read_classical_blocks(path):
    with h5py.File(path, "r") as f:
        g = f["readout"]
        cell = lambda n: ["".join(chr(int(c)) for c in np.array(f[r]).ravel()) for r in np.array(g[n]).ravel()]
        dom, roi = cell("domains"), cell("roiNames")
        locRate = np.array(g["locRate"]).astype(float).transpose()   # [nDom,nRung,maxBlk,nROI]
        bgLoc = np.array(g["bgLoc"]).astype(float).transpose()        # [nDom,nBg,nROI]
    d = dom.index("PI"); nRung = locRate.shape[1]

    def blocks(rr_name):
        rr = roi.index(rr_name)
        rung = [locRate[d, i, :, rr][np.isfinite(locRate[d, i, :, rr])] for i in range(nRung)]
        bg = bgLoc[d, :, rr][np.isfinite(bgLoc[d, :, rr])]
        return rung, bg
    return blocks


def read_learned_blocks(npy_path, method, cfg):
    out = np.load(npy_path, allow_pickle=True).item()
    rung_tags = DS.rung_blocks(cfg); bg_tags = list(cfg["bgFlow"])

    def blocks(rr_name):
        pb = out["results"]["PI"][rr_name][method]["locRate"]
        rung = [np.array([pb[t] for t in tags], float) for tags in rung_tags]
        bg = np.array([pb[t] for t in bg_tags], float)
        return rung, bg
    return blocks


def loc_sub_from_blocks(rung, bg):
    lm = np.array([np.mean(r) if len(r) else np.nan for r in rung])
    return np.maximum(lm - np.mean(bg), 0.0)


def win_slope(conc, idxs, ls):
    idxs = np.array(idxs); x = np.log10(conc[idxs]); y = ls[idxs]
    m = np.isfinite(x) & np.isfinite(y) & (y > 0)
    if m.sum() < 2:
        return np.nan
    return np.polyfit(x[m], np.log10(y[m]), 1)[0]


def breakpoint(conc, idxs, ls):
    idxs = np.array(idxs); x = np.log10(conc[idxs]); y = ls[idxs]
    m = np.isfinite(x) & (y > 0); x = x[m]; y = np.log10(y[m])
    o = np.argsort(x); x, y = x[o], y[o]
    if x.size < 6:
        return (np.nan, np.nan, np.nan)
    grid = np.linspace(x[2], x[-3], 80); best = None
    for psi in grid:
        X = np.column_stack([np.ones_like(x), x, np.maximum(x - psi, 0)])
        b, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
        sse = np.sum((X @ b - y) ** 2)
        if best is None or sse < best[0]:
            best = (sse, psi, b)
    _, psi, b = best
    return (10 ** psi, b[1], b[1] + b[2])


def sliding_windows(conc, exclude, w):
    valid = sorted([i for i in range(len(conc)) if i not in exclude], key=lambda i: conc[i])
    return [valid[k:k + w] for k in range(len(valid) - w + 1)]


def linear_range_end(conc, swins, ls):
    """Highest-concentration sliding window whose slope is still >= 1; return its top concentration.
    NaN if no window reaches linear."""
    for win in reversed(swins):                 # highest-concentration window first
        s = win_slope(conc, win, ls)
        if np.isfinite(s) and s >= 1.0:
            return float(np.max(conc[win]))
    return np.nan


def _ci(a):
    a = np.asarray(a, float); a = a[np.isfinite(a)]
    return (np.percentile(a, 2.5), np.percentile(a, 97.5)) if a.size else (np.nan, np.nan)


def analyze(conc, windows, exclude, slide_w, blocks_fn, roi_name):
    rung, bg = blocks_fn(roi_name)
    ls = loc_sub_from_blocks(rung, bg)
    swins = sliding_windows(conc, exclude, slide_w)
    all_idx = windows["all"]
    pt = {w: win_slope(conc, idxs, ls) for w, idxs in windows.items()}
    bp = breakpoint(conc, all_idx, ls)
    lre = linear_range_end(conc, swins, ls)
    rng = np.random.default_rng(SEED)
    boot = {w: [] for w in windows}; bpc = []; bsb = []; bsa = []; lreb = []
    for _ in range(NBOOT):
        rb = [(r[rng.integers(0, len(r), len(r))] if len(r) else r) for r in rung]
        bb = bg[rng.integers(0, len(bg), len(bg))]
        lsb = loc_sub_from_blocks(rb, bb)
        for w, idxs in windows.items():
            boot[w].append(win_slope(conc, idxs, lsb))
        c, sb, sa = breakpoint(conc, all_idx, lsb); bpc.append(c); bsb.append(sb); bsa.append(sa)
        lreb.append(linear_range_end(conc, swins, lsb))
    return dict(loc_sub=ls, pt=pt, ci={w: _ci(boot[w]) for w in windows},
                bp=bp, bp_ci=_ci(bpc), sb_ci=_ci(bsb), sa_ci=_ci(bsa),
                lre=lre, lre_ci=_ci(lreb))


def main():
    Z = OUTPUT_ROOT                                     # reanalysis_out (readout_<ds>.mat)
    SC = sys.argv[1] if len(sys.argv) > 1 else SCRATCH_DIR  # holds locrate_<ds>.npy
    MK = {"Classical": "classical", "Deep-ULM": "deepulm", "LISTA": "lista"}
    stats = {}
    for ds in ("jun23", "apr17"):
        dc = DSC[ds]; cfg = DS.get_dataset(ds); conc = np.asarray(cfg["conc"], float)
        hi_win = "U1-U5" if ds == "jun23" else "C5-C8"
        classical = read_classical_blocks(os.path.join(Z, ds, "readout_%s.mat" % ds))
        methods = {"Classical": classical,
                   "Deep-ULM": read_learned_blocks(os.path.join(SC, "locrate_%s.npy" % ds), "deepulm", cfg),
                   "LISTA": read_learned_blocks(os.path.join(SC, "locrate_%s.npy" % ds), "lista", cfg)}
        print("\n" + "#" * 112)
        print("DATASET: %s    (%s)" % (ds, "CEILING TEST" if ds == "jun23" else "NEGATIVE CONTROL: below the ceiling, expect no clean breakpoint"))
        stats[ds] = {}
        for tube in dc["tubes"]:
            print("=" * 112)
            print("TUBE: %s" % tube)
            print("%-10s | %-24s | %-6s | %-6s | %-18s | %-22s" %
                  ("method", "breakpoint C [95% CI]", "s_bel", "s_abv", "%s slope [95%% CI]" % hi_win,
                   "linear-range-end C [95% CI]"))
            print("-" * 112)
            stats[ds][tube] = {}
            for mname, fn in methods.items():
                r = analyze(conc, dc["windows"], dc["exclude"], dc["slide_w"], fn, tube)
                bp, bpc = r["bp"], r["bp_ci"]; h, hci = r["pt"][hi_win], r["ci"][hi_win]
                lre, lci = r["lre"], r["lre_ci"]
                print("%-10s | %8.2e [%.1e,%.1e] | %5.2f | %5.2f | %5.2f [%5.2f,%5.2f] | %8.2e [%.1e,%.1e]" % (
                    mname, bp[0], bpc[0], bpc[1], bp[1], bp[2], h, hci[0], hci[1], lre, lci[0], lci[1]))
                stats[ds][tube][MK[mname]] = dict(
                    breakpoint=bp[0], breakpoint_ci=bpc, slope_below=bp[1], slope_above=bp[2],
                    u1u5=h, u1u5_ci=hci, lre=lre, lre_ci=lci,
                    windows={w: (r["pt"][w], r["ci"][w]) for w in dc["windows"]})
            print("  per-window slopes [95% CI]:")
            for mname, fn in methods.items():
                r = analyze(conc, dc["windows"], dc["exclude"], dc["slide_w"], fn, tube)
                s = "  ".join("%s %.2f[%.2f,%.2f]" % (w, r["pt"][w], r["ci"][w][0], r["ci"][w][1]) for w in dc["windows"])
                print("    %-10s %s" % (mname, s))
    np.save(os.path.join(SC, "ceiling_stats.npy"), stats, allow_pickle=True)
    print("\n[ceiling] stats saved -> %s" % os.path.join(SC, "ceiling_stats.npy"))


if __name__ == "__main__":
    main()
