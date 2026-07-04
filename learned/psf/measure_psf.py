#!/usr/bin/env python
# SPDX-License-Identifier: MIT
"""
Phase 3 PSF measurement: rough anisotropic-Gaussian PSF from isolated single bubbles.

  >>> STOPGAP PSF <<<  See phase3/UPGRADE_LATER.md.
  Decision (Eli, 2026-06-18): sanity-check + reuse. Measure a ROUGH anisotropic Gaussian from
  apr17 single bubbles (prefer C3-C6; check whether C1/C2 are cleaner), VALIDATE against the
  analytical model psf.mat (lat FWHM 0.2447 mm, ax FWHM 0.0694 mm), train Phase 3 on the
  analytical / model-refined PSF NOW. UPGRADE to a measured LINEAR PSF from the followup later.

Method, per block:
  1. read N frames of PI-sum compound IQ (HDF5, complex) from the Z: cache
  2. SVD clutter filter (remove top-k=8 singular vectors), envelope=|.|
  3. detect bright, ISOLATED local maxima INSIDE THE TUBE ROI (where bubbles actually are)
  4. lateral measurement: 1D Gaussian fit to the lateral profile through the peak (well-sampled,
     FWHM~7px) -> reliable lateral FWHM. quality gate = 1D lateral r2.
  5. axial: 1D fit to the axial profile (UNDERSAMPLED, FWHM~2px = sampling floor) -> reported
     as secondary / lower bound only, NOT a trustworthy measurement on this grid.

Axis convention: HDF5 IQ_compound is (nFrames, nX=193, nZ=404); each frame transposed to
[nZ=404 axial, nX=193 lateral]. dx=dz=0.0346875 mm. Tube ROI x[-7.5,-5.25], z[16,22] mm.

Usage:
  python measure_psf.py --blocks C4b1 --frames 256
  python measure_psf.py --blocks C1b1 C2b1 C3b1 C4b1 C5b1 C6b1 --frames 512
"""
import os, sys, time, argparse, json
import numpy as np
import h5py
from scipy.sparse.linalg import svds
from scipy.optimize import curve_fit
from scipy.ndimage import maximum_filter

BASE = r"Z:\Eli\4-17-FINAL\IQ_cache"
DX_MM = DZ_MM = 0.0346875
FWHM_K = 2.0 * np.sqrt(2.0 * np.log(2.0))           # 2.3548
MODEL_FWHM_AX_MM, MODEL_FWHM_LAT_MM = 0.069375, 0.24470898
# grid extents (from batch_config.mat): xGrid[-9.095,-2.435] (193), zGrid[13.951,27.93] (404)
XG0, ZG0 = -9.095, 13.951
TUBE_X_MM, TUBE_Z_MM = (-7.5, -5.25), (16.0, 22.0)


def mm_to_ix(mm, g0):
    return int(round((mm - g0) / DX_MM))


TUBE = dict(x0=mm_to_ix(TUBE_X_MM[0], XG0), x1=mm_to_ix(TUBE_X_MM[1], XG0),
            z0=mm_to_ix(TUBE_Z_MM[0], ZG0), z1=mm_to_ix(TUBE_Z_MM[1], ZG0))


def resolve_path(tag):
    needle = tag + "-" if tag[0] == "C" else tag
    hits = [f for f in os.listdir(BASE) if (needle in f) and f.endswith("_IQ.mat")]
    if not hits:
        raise FileNotFoundError(f"no _IQ.mat for {tag!r}")
    if len(hits) > 1:
        raise RuntimeError(f"ambiguous {tag!r}: {hits}")
    return os.path.join(BASE, hits[0])


def read_block(tag, nframes):
    with h5py.File(resolve_path(tag), "r") as f:
        d = f["IQ_compound"]
        F = min(nframes, d.shape[0])
        sl = d[:F]
    c = sl["real"].astype(np.float32) + 1j * sl["imag"].astype(np.float32)   # [F,193,404]
    return np.ascontiguousarray(np.transpose(c, (0, 2, 1)))                   # [F,404,193]


def svd_clutter_filter(stack, k=8):
    F, nZ, nX = stack.shape
    M = stack.reshape(F, nZ * nX).T
    U, s, Vh = svds(M.astype(np.complex64), k=k, which="LM")
    filt = M - (U * s) @ Vh
    return np.abs(filt).T.reshape(F, nZ, nX).astype(np.float32)


def _g1d(x, A, mu, sig, off):
    return A * np.exp(-((x - mu) ** 2) / (2 * sig ** 2)) + off


def fit_1d(profile, c, half):
    """Fit 1D Gaussian to profile[c-half:c+half+1]; index space relative to c. Returns sig,r2,off."""
    lo, hi = c - half, c + half + 1
    y = profile[lo:hi].astype(np.float64)
    x = np.arange(len(y), dtype=np.float64)
    off0 = np.median(y)
    A0 = y.max() - off0
    try:
        popt, _ = curve_fit(_g1d, x, y, p0=[A0, half, 2.0, off0],
                            bounds=([0, half - 2, 0.3, -np.inf], [np.inf, half + 2, 12, np.inf]),
                            maxfev=3000)
    except Exception:
        return None
    resid = y - _g1d(x, *popt)
    ss_res = float(np.sum(resid ** 2)); ss_tot = float(np.sum((y - y.mean()) ** 2)) + 1e-12
    return dict(sig=float(popt[2]), r2=1.0 - ss_res / ss_tot, A=float(popt[0]), off=float(popt[3]))


def detect_isolated_roi(frame, thr, iso_px, wlat, wax):
    mx = maximum_filter(frame, size=3)
    peaks = (frame == mx) & (frame > thr)
    zz, xx = np.where(peaks)
    keep = (zz >= TUBE["z0"]) & (zz <= TUBE["z1"]) & (xx >= TUBE["x0"]) & (xx <= TUBE["x1"]) \
        & (zz >= wax + 1) & (zz < frame.shape[0] - wax - 1) \
        & (xx >= wlat + 1) & (xx < frame.shape[1] - wlat - 1)
    pts = np.column_stack([zz[keep], xx[keep]])
    if len(pts) == 0:
        return pts
    iso = []
    for i, (z, x) in enumerate(pts):
        d2 = ((pts[:, 0] - z) ** 2 + (pts[:, 1] - x) ** 2).astype(float); d2[i] = np.inf
        if np.sqrt(d2.min()) > iso_px:
            iso.append((z, x))
    return np.array(iso, dtype=int) if iso else np.empty((0, 2), int)


def process_block(tag, nframes, pct=99.9, iso_px=14, wlat=12, wax=6, r2_min=0.90, k=8):
    t0 = time.time(); stack = read_block(tag, nframes); t_read = time.time() - t0
    t1 = time.time(); env = svd_clutter_filter(stack, k=k); t_svd = time.time() - t1
    thr = float(np.percentile(env, pct))
    n_cand = 0; recs = []
    for fr in range(env.shape[0]):
        f = env[fr]
        pts = detect_isolated_roi(f, thr, iso_px, wlat, wax)
        n_cand += len(pts)
        for z, x in pts:
            lat = fit_1d(f[z, :], x, wlat)              # lateral profile through peak row
            ax = fit_1d(f[:, x], z, wax)                # axial profile through peak col
            if lat is None:
                continue
            recs.append(dict(frame=int(fr), z=int(z), x=int(x), amp=float(f[z, x]),
                             s_lat=lat["sig"], r2_lat=lat["r2"],
                             s_ax=(ax["sig"] if ax else np.nan),
                             r2_ax=(ax["r2"] if ax else np.nan)))
    good = [r for r in recs if r["r2_lat"] >= r2_min and 1.0 < r["s_lat"] < 8.0]
    sl = np.array([r["s_lat"] for r in good]); sa = np.array([r["s_ax"] for r in good])
    def med_fwhm(arr, dd): return float(np.nanmedian(arr) * FWHM_K * dd) if len(arr) else None
    summ = dict(tag=tag, nframes=int(env.shape[0]), thr=thr, n_candidates=int(n_cand),
                n_good=len(good),
                fwhm_lat_mm=med_fwhm(sl, DX_MM),
                fwhm_lat_iqr=float((np.percentile(sl, 75) - np.percentile(sl, 25)) * FWHM_K * DX_MM) if len(sl) else None,
                fwhm_ax_mm=med_fwhm(sa, DZ_MM),
                med_r2_lat=float(np.median([r["r2_lat"] for r in good])) if good else None,
                t_read=t_read, t_svd=t_svd)
    return summ, good


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--blocks", nargs="+", required=True)
    ap.add_argument("--frames", type=int, default=512)
    ap.add_argument("--pct", type=float, default=99.9)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    print(f"[psf] tube ROI px: z[{TUBE['z0']},{TUBE['z1']}] x[{TUBE['x0']},{TUBE['x1']}]")
    print(f"[psf] model FWHM_lat={MODEL_FWHM_LAT_MM:.4f}  FWHM_ax={MODEL_FWHM_AX_MM:.4f} mm")
    print(f"[psf] {'tag':6} {'nfr':>4} {'cand':>5} {'good':>5} {'FWHMlat':>8} {'+-IQR':>6} "
          f"{'FWHMax*':>7} {'r2lat':>6}  read/svd")
    summ_all, good_all = [], {}
    for tag in args.blocks:
        s, g = process_block(tag, args.frames, pct=args.pct)
        summ_all.append(s); good_all[tag] = g
        fl = f"{s['fwhm_lat_mm']:.4f}" if s['fwhm_lat_mm'] else "  --  "
        iqr = f"{s['fwhm_lat_iqr']:.3f}" if s['fwhm_lat_iqr'] else " -- "
        fa = f"{s['fwhm_ax_mm']:.4f}" if s['fwhm_ax_mm'] else "  --  "
        r2 = f"{s['med_r2_lat']:.3f}" if s['med_r2_lat'] else " -- "
        print(f"[psf] {tag:6} {s['nframes']:>4} {s['n_candidates']:>5} {s['n_good']:>5} "
              f"{fl:>8} {iqr:>6} {fa:>7} {r2:>6}  {s['t_read']:.0f}/{s['t_svd']:.0f}")
    # pooled lateral
    pooled = np.concatenate([[r["s_lat"] for r in good_all[t]] for t in good_all if good_all[t]]) \
        if any(good_all.values()) else np.array([])
    if len(pooled):
        flat = np.median(pooled) * FWHM_K * DX_MM
        print(f"\n[psf] POOLED lateral n={len(pooled)}  FWHM_lat={flat:.4f} mm  "
              f"(model {MODEL_FWHM_LAT_MM:.4f}, ratio {flat/MODEL_FWHM_LAT_MM:.2f})")
        print("[psf] axial FWHM is at the ~2px sampling floor -> NOT a reliable measurement; "
              "model uses lambda (0.0694 mm). See UPGRADE_LATER.md.")
    out = args.out or os.path.join(os.path.dirname(__file__), "..", "artifacts",
                                   "psf_measure_" + "_".join(args.blocks) + ".npz")
    np.savez(out, summary=json.dumps(summ_all),
             **{f"{t}_slat": np.array([r["s_lat"] for r in good_all[t]]) for t in good_all},
             **{f"{t}_sax": np.array([r["s_ax"] for r in good_all[t]]) for t in good_all},
             **{f"{t}_r2": np.array([r["r2_lat"] for r in good_all[t]]) for t in good_all})
    print(f"[psf] saved -> {os.path.abspath(out)}")


if __name__ == "__main__":
    main()
