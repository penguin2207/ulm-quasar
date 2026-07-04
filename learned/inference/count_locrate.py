#!/usr/bin/env python
# SPDX-License-Identifier: MIT
"""
Count: in-tube LOCALIZATION RATE from external (Python-net) pixel localizations.
There is NO tracking and NO QC. The count is, per (dataset, domain, method):

  per block, per ROI:
    locRate = (# detections whose pixel [row,col] land in the ROI mask, over ALL frames) / nFrames
  per rung:
    loc_mean = mean(locRate over that rung's blocks)            (variable blocks-per-rung)
    pedLoc   = mean(locRate over the bgFlow blocks)
    loc_sub  = max(loc_mean - pedLoc, 0)                        (aggregate-then-subtract, floor 0)
  slope:
    lf fit (Section 6c): log10-log10 polyfit, 400-boot percentile CI, NO weights, NO Poisson.

ROI masks are rasterized per block from roi_polys_<ds>.mat (var 'roi') via matplotlib.path.Path
(the inpolygon analogue, Section 4b). roi.names is the FROZEN 1x5 order
{full, combinedTube, tubeL, tubeR, background}; roi.poly{1}=[] ('full') is the whole-FOV sentinel.

This module is import-safe: the pure functions (rasterize_roi, count_in_roi, loc_rate, aggregate
loc_sub, lf_fit, beta keying) are unit-tested without any real caches. The disk-driven driver
(run_dataset) needs the STEP-2 outputs (per-block loc .mats + roi_polys_<ds>.mat) to run.
"""
import os, sys, argparse, warnings
import numpy as np
from matplotlib.path import Path

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import datasets as DS                                            # noqa: E402


# ============================================================ pure / testable core
def rasterize_roi(poly, xGrid, zGrid):
    """Rasterize one ROI polygon onto a block grid -> [nZ x nX] boolean mask (Section 4b).

    poly : [N x 2] physical vertices [x_mm, z_mm], OR None/empty for the 'full' sentinel
           (whole-FOV: mask all true).
    Uses matplotlib.path.Path(poly).contains_points over meshgrid(xGrid, zGrid), the spec's
    inpolygon analogue (inpolygon(XX, ZZ, poly(:,1), poly(:,2))).
    """
    nZ, nX = len(zGrid), len(xGrid)
    if poly is None or len(np.atleast_2d(poly)) == 0 or np.size(poly) == 0:
        return np.ones((nZ, nX), dtype=bool)                    # 'full' sentinel
    poly = np.asarray(poly, dtype=np.float64).reshape(-1, 2)
    XX, ZZ = np.meshgrid(xGrid, zGrid)                          # both [nZ x nX]
    pts = np.column_stack([XX.ravel(), ZZ.ravel()])            # (x_mm, z_mm)
    mask = Path(poly).contains_points(pts).reshape(nZ, nX)
    return mask


def count_in_roi(rows, cols, mask):
    """# detections whose NATIVE pixel [row, col] (0-based) land inside `mask` [nZ x nX]."""
    rows = np.asarray(rows).astype(np.int64).ravel()
    cols = np.asarray(cols).astype(np.int64).ravel()
    if rows.size == 0:
        return 0
    nZ, nX = mask.shape
    ok = (rows >= 0) & (rows < nZ) & (cols >= 0) & (cols < nX)
    if not ok.any():
        return 0
    return int(mask[rows[ok], cols[ok]].sum())


def loc_rate(rows, cols, mask, nFrames):
    """In-ROI localization rate = (in-ROI detections over ALL frames) / nFrames (Section 2 step 5)."""
    if nFrames is None or nFrames <= 0:
        return np.nan
    return count_in_roi(rows, cols, mask) / float(nFrames)


def aggregate_loc_sub(block_rates_by_rung, ped):
    """loc_mean per rung (mean over that rung's blocks, omitnan) and
    loc_sub = max(loc_mean - ped, 0) (Section 7, aggregate-then-subtract, floor at 0).

    block_rates_by_rung : list (len nRung) of array-likes of per-block rates (variable length).
    ped                 : scalar pedestal (mean over bgFlow blocks).
    Returns (loc_mean [nRung], loc_sub [nRung])."""
    nRung = len(block_rates_by_rung)
    loc_mean = np.full(nRung, np.nan)
    for r in range(nRung):
        vals = np.asarray(block_rates_by_rung[r], dtype=float).ravel()
        vals = vals[np.isfinite(vals)]
        if vals.size:
            loc_mean[r] = vals.mean()
    ped = 0.0 if (ped is None or not np.isfinite(ped)) else float(ped)
    loc_sub = np.maximum(loc_mean - ped, 0.0)
    return loc_mean, loc_sub


def pedestal(bg_rates):
    """Bg pedestal = mean over bgFlow blocks' rates (omitnan); empty -> 0 (Section 7)."""
    vals = np.asarray(bg_rates, dtype=float).ravel()
    vals = vals[np.isfinite(vals)]
    return float(vals.mean()) if vals.size else 0.0


def lf_fit(x, y, n_boot=400, seed=0):
    """`lf` local fit (RESPONSE.md Section 6c): mask isfinite & x>0 & y>0; need n>=3 else NaN.
    X=log10(x), Y=log10(y), beta = polyfit(X,Y,1)[0]; R2 = 1-SSE/SST; bootstrap n_boot=400
    percentile CI [.025,.975]. NO weighting, NO Poisson se.

    NOTE: RESPONSE.md 6c does not state a bootstrap RNG seed; we fix one for reproducibility
    (the point estimate `beta` and `R2` are seed-independent; only `ci` depends on it)."""
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    m = np.isfinite(x) & np.isfinite(y) & (x > 0) & (y > 0)
    x, y = x[m], y[m]
    n = int(x.size)
    if n < 3:
        return dict(n=n, beta=float("nan"), ci=[float("nan"), float("nan")], R2=float("nan"))
    X, Y = np.log10(x), np.log10(y)
    p = np.polyfit(X, Y, 1)
    beta = float(p[0])
    yhat = np.polyval(p, X)
    sse = float(np.sum((Y - yhat) ** 2))
    sst = float(np.sum((Y - Y.mean()) ** 2))
    R2 = (1.0 - sse / sst) if sst > 0 else float("nan")
    rng = np.random.default_rng(seed)
    bs = np.empty(n_boot)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")                       # collinear resamples -> harmless RankWarning
        for b in range(n_boot):
            idx = rng.integers(0, n, n)                        # MATLAB randi(n,n,1) analogue
            bs[b] = np.polyfit(X[idx], Y[idx], 1)[0]
    ci = np.quantile(bs, [0.025, 0.975])
    return dict(n=n, beta=beta, ci=[float(ci[0]), float(ci[1])], R2=float(R2))


def make_beta_key(domain, roi, metric, window):
    """betas.lookup key = matlab.lang.makeValidName(sprintf('%s_%s_%s_%s',...)) analogue
    (RESPONSE.md Section 6b). For these keys the only transform is '-' -> '_' in the window
    label (e.g. C3-C7 -> C3_C7); we also coerce any other non [A-Za-z0-9_] to '_'."""
    raw = "%s_%s_%s_%s" % (domain, roi, metric, window)
    out = "".join(c if (c.isalnum() or c == "_") else "_" for c in raw)
    if out and out[0].isdigit():
        out = "x" + out
    return out


# ============================================================ disk loaders (need STEP-2 data)
def load_roi(roi_file):
    """Load roi_polys_<ds>.mat (var 'roi') -> (names[list], polys[list]); poly[0] ('full') is
    None (whole-FOV sentinel). RESPONSE.md Section 4b."""
    from scipy.io import loadmat
    roi = loadmat(roi_file, squeeze_me=True, struct_as_record=False)["roi"]
    names = [str(s) for s in np.atleast_1d(roi.names)]
    raw_polys = np.atleast_1d(roi.poly)
    polys = []
    for rr in range(len(names)):
        p = raw_polys[rr]
        p = None if (p is None or np.size(p) == 0) else np.asarray(p, float).reshape(-1, 2)
        polys.append(p)
    return names, polys


def load_block_locs(loc_path):
    """Load one infer output .mat -> (rows, cols, amps, frames, grid_dict, nFrames)."""
    from scipy.io import loadmat
    m = loadmat(loc_path, squeeze_me=True, struct_as_record=False)
    L = np.atleast_2d(m["localizations"])
    if L.size == 0:
        L = L.reshape(0, 4)
    g = m["g"]
    grid = dict(xGrid=np.atleast_1d(g.xGrid).astype(float).ravel(),
                zGrid=np.atleast_1d(g.zGrid).astype(float).ravel(),
                nZ=int(g.nZ), nX=int(g.nX))
    nFrames = int(np.atleast_1d(m["nFrames"]).ravel()[0])
    rows = L[:, 0]; cols = L[:, 1]; amps = L[:, 2]; frames = L[:, 3]
    return rows, cols, amps, frames, grid, nFrames


def _loc_path(loc_dir, method, tag, domain):
    return os.path.join(loc_dir, "%s_%s_%s_locs.mat" % (method, tag, domain))


# ============================================================ driver
def run_dataset(dataset, loc_dir, out_path=None, domains=None, methods=("deconv", "deepulm", "lista"),
                roi_file=None, verbose=True):
    """Build results[domain][roi][method] = {locRate per block, loc_mean, pedLoc, loc_sub, fit}.

    Reads per-block infer outputs '<method>_<tag>_<domain>_locs.mat' from `loc_dir`, rasterizes the
    ROI polygons per block grid, counts in-ROI detections / nFrames, aggregates with the bgFlow
    pedestal, and lf-fits loc_sub vs conc on the headline + full + all configured windows.
    """
    cfg = DS.get_dataset(dataset) if isinstance(dataset, str) else dataset
    if domains is None:
        domains = list(cfg["domains"])
    if roi_file is None:
        roi_file = cfg["roiFile"]
    roi_names, roi_polys = load_roi(roi_file)

    rung_tag_lists = DS.rung_blocks(cfg)                        # variable blocks per rung
    bg_tags = list(cfg["bgFlow"])
    conc = np.asarray(cfg["conc"], float)

    results = {}
    for domain in domains:
        results[domain] = {}
        for method in methods:
            # per-block rates: roi -> list-per-rung of block rates ; roi -> list of bg rates
            rung_rates = {rr: [[] for _ in range(cfg["nRung"])] for rr in roi_names}
            bg_rates = {rr: [] for rr in roi_names}
            per_block = {rr: {} for rr in roi_names}           # roi -> {tag: rate}
            for ri, tags in enumerate(rung_tag_lists):
                for tag in tags:
                    rr_rate = _one_block_rates(loc_dir, method, tag, domain, roi_names, roi_polys, verbose)
                    for rr in roi_names:
                        rung_rates[rr][ri].append(rr_rate[rr])
                        per_block[rr][tag] = rr_rate[rr]
            for tag in bg_tags:
                rr_rate = _one_block_rates(loc_dir, method, tag, domain, roi_names, roi_polys, verbose)
                for rr in roi_names:
                    bg_rates[rr].append(rr_rate[rr])
                    per_block[rr][tag] = rr_rate[rr]
            # ---- aggregate + fit per ROI ----
            for rr in roi_names:
                ped = pedestal(bg_rates[rr])
                loc_mean, loc_sub = aggregate_loc_sub(rung_rates[rr], ped)
                fits = {}
                for window in cfg["fitWindows"]:
                    lo, hi = window
                    wl = DS.window_label(cfg, window)
                    sel = slice(lo - 1, hi)
                    fit_sub = lf_fit(conc[sel], loc_sub[sel])
                    fit_raw = lf_fit(conc[sel], loc_mean[sel])
                    fits[wl] = dict(window=window, is_headline=(window == cfg["headlineWindow"]),
                                    beta_sub=fit_sub["beta"], ci_sub=fit_sub["ci"], R2_sub=fit_sub["R2"],
                                    beta_raw=fit_raw["beta"], ci_raw=fit_raw["ci"], R2_raw=fit_raw["R2"],
                                    n=fit_sub["n"], beta_key=make_beta_key(domain, rr, "locRate", wl))
                results[domain].setdefault(rr, {})[method] = dict(
                    locRate=per_block[rr], loc_mean=loc_mean, pedLoc=ped, loc_sub=loc_sub,
                    conc=conc, fits=fits)

    out = dict(dataset=cfg["name"], domains=list(domains), roiNames=roi_names,
               rungLabels=list(cfg["rungLabels"]), conc=conc,
               headlineDomain=cfg["headlineDomain"], headlineROI=cfg["headlineROI"],
               headlineWindow=cfg["headlineWindow"], results=results)
    if out_path:
        np.save(out_path, out, allow_pickle=True)
        if verbose:
            print("[count] saved -> %s" % out_path)
    return out


def _one_block_rates(loc_dir, method, tag, domain, roi_names, roi_polys, verbose):
    """Per-ROI loc rate for one (method, tag, domain) block. Missing file -> all NaN."""
    path = _loc_path(loc_dir, method, tag, domain)
    if not os.path.exists(path):
        if verbose:
            print("    [%s %s %s] MISSING -> NaN" % (method, tag, domain))
        return {rr: np.nan for rr in roi_names}
    rows, cols, amps, frames, grid, nFrames = load_block_locs(path)
    out = {}
    for rr_name, poly in zip(roi_names, roi_polys):
        mask = rasterize_roi(poly, grid["xGrid"], grid["zGrid"])
        out[rr_name] = loc_rate(rows, cols, mask, nFrames)
    return out


# ============================================================ anchor readers (Section 6)
def read_anchor_readout(readout_file):
    """Read readout_<ds>.mat (var readout): loc_sub [nDom x nRung x nROI], domains, roiNames,
    rungConc, thrFixed (Section 6a). Handles BOTH -v7 (scipy) and -v7.3/HDF5 (h5py): the real
    RUNNER writes readout as -v7.3, which scipy.loadmat cannot read. Returns a dict of the
    count-overlay-relevant fields."""
    try:
        from scipy.io import loadmat
        r = loadmat(readout_file, squeeze_me=True, struct_as_record=False)["readout"]
    except NotImplementedError:
        return _read_anchor_readout_h5(readout_file)
    out = dict(domains=[str(s) for s in np.atleast_1d(r.domains)],
               roiNames=[str(s) for s in np.atleast_1d(r.roiNames)],
               rungConc=np.atleast_1d(r.rungConc).astype(float).ravel(),
               loc_sub=np.asarray(r.loc_sub, float))
    for f in ("rungLabels", "thrFixed", "loc_mean", "pedLoc", "locRate"):
        if hasattr(r, f):
            v = getattr(r, f)
            out[f] = [str(s) for s in np.atleast_1d(v)] if f == "rungLabels" else np.asarray(v, float)
    return out


def _read_anchor_readout_h5(readout_file):
    """v7.3/HDF5 fallback for read_anchor_readout: dereference MATLAB cellstr and reverse array
    axes (HDF5 row-major vs MATLAB column-major). Validated against the real readout_apr17.mat:
    loc_sub h5 (nROI,nRung,nDom) -> transpose -> (nDom,nRung,nROI); thrFixed -> [418.97,423.40,294.61]."""
    import h5py
    with h5py.File(readout_file, "r") as f:
        g = f["readout"]

        def cellstr(name):
            return ["".join(chr(int(c)) for c in np.array(f[ref]).ravel())
                    for ref in np.array(g[name]).ravel()]

        out = dict(domains=cellstr("domains"),
                   roiNames=cellstr("roiNames"),
                   rungConc=np.array(g["rungConc"]).astype(float).ravel(),
                   loc_sub=np.array(g["loc_sub"]).astype(float).transpose())
        if "thrFixed" in g:
            out["thrFixed"] = np.array(g["thrFixed"]).astype(float).ravel()
        for fld in ("loc_mean", "pedLoc"):
            if fld in g:
                out[fld] = np.array(g[fld]).astype(float).transpose()
        if "rungLabels" in g:
            out["rungLabels"] = cellstr("rungLabels")
    return out


def read_anchor_betas(betas_file):
    """Read betas_<ds>.mat (var betas): betas.lookup.(key) structs (Section 6b). Returns
    {key: {beta_sub, ci_sub, R2_sub, beta_raw, ci_raw, R2_raw, domain, roi, metric, window, n}}."""
    from scipy.io import loadmat
    betas = loadmat(betas_file, squeeze_me=True, struct_as_record=False)["betas"]
    lookup = betas.lookup
    out = {}
    for key in (lookup._fieldnames if hasattr(lookup, "_fieldnames") else []):
        s = getattr(lookup, key)
        out[key] = {f: _coerce(getattr(s, f)) for f in s._fieldnames}
    return out


def anchor_count_axis(readout, domain="PI", roi="combinedTube"):
    """Anchor loc_sub vs conc for a (domain, roi): returns (x=rungConc, y=loc_sub[d,:,rr])."""
    d = readout["domains"].index(domain)
    rr = readout["roiNames"].index(roi)
    return readout["rungConc"], np.asarray(readout["loc_sub"])[d, :, rr]


def _coerce(v):
    if hasattr(v, "ravel") and np.size(v) > 1:
        return np.asarray(v).ravel().tolist()
    try:
        return float(v)
    except (TypeError, ValueError):
        return str(v)


# ============================================================ CLI
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=list(DS.DATASETS))
    ap.add_argument("--loc_dir", required=True, help="dir with <method>_<tag>_<domain>_locs.mat")
    ap.add_argument("--out", default=None, help="output .npy (pickled results dict)")
    ap.add_argument("--domains", default=None, help="comma list (default: all 3)")
    ap.add_argument("--methods", default="deconv,deepulm,lista")
    ap.add_argument("--roi_file", default=None, help="roi_polys_<ds>.mat (default: from config)")
    args = ap.parse_args()
    domains = args.domains.split(",") if args.domains else None
    methods = tuple(args.methods.split(","))
    run_dataset(args.dataset, args.loc_dir, out_path=args.out, domains=domains,
                methods=methods, roi_file=args.roi_file)


if __name__ == "__main__":
    main()
