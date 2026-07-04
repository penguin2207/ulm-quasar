#!/usr/bin/env python
# SPDX-License-Identifier: MIT
"""
Phase 3 classical baseline: sparse deconvolution by FISTA (van Sloun 2021 Sec II.D.2).

  minimize  0.5 || y - A x ||^2 + lambda || x ||_1 ,  x >= 0
  A x = downsample_s( conv(x_SR, PSF_SR) )         (shifted-PSF dictionary, x8 up-sampling)
  lambda = 0.01, x8 grid up-sampling (paper spec).

No training. Stands up the inference + peak-extraction + counting scaffold the learned methods
reuse. Implemented in torch (device-agnostic: CPU now, GPU on the 3070 once a CUDA build is
installed). Uses the STOPGAP PSF (phase3/artifacts/psf_training.npz) -- see UPGRADE_LATER.md.

The forward/adjoint pair is exact-adjoint (verified by power-iteration Lipschitz); avg-pool is
the native imaging integral, its adjoint is nearest-upsample / s^2.
"""
import os, sys, time, json, argparse
import numpy as np
import torch
import torch.nn.functional as F

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.normpath(os.path.join(HERE, "..")))
from datagen.synth import load_psf, make_sample, SynthConfig          # noqa: E402
from common.peaks import extract_peaks_sr, match_localizations         # noqa: E402


def render_psf_sr(sx_px, sz_px, sr, device="cpu", nsig=4.0):
    sxs, szs = sx_px * sr, sz_px * sr
    hz, hx = int(np.ceil(nsig * szs)), int(np.ceil(nsig * sxs))
    zz, xx = np.mgrid[-hz:hz + 1, -hx:hx + 1]
    k = np.exp(-(zz ** 2) / (2 * szs ** 2) - (xx ** 2) / (2 * sxs ** 2))
    k /= k.sum()
    return torch.tensor(k, dtype=torch.float32, device=device)


def _conv(x, k):    # true convolution, 'same'
    return F.conv2d(x[None, None], k.flip(0, 1)[None, None],
                    padding=(k.shape[0] // 2, k.shape[1] // 2))[0, 0]


def _corr(x, k):    # correlation = adjoint of convolution
    return F.conv2d(x[None, None], k[None, None],
                    padding=(k.shape[0] // 2, k.shape[1] // 2))[0, 0]


def forward_A(x_sr, k, sr):
    return F.avg_pool2d(_conv(x_sr, k)[None, None], sr)[0, 0]


def adjoint_At(r, k, sr):
    up = F.interpolate(r[None, None], scale_factor=sr, mode="nearest")[0, 0] / (sr * sr)
    return _corr(up, k)


def lipschitz(k, sr, shape_sr, device, iters=12):
    v = torch.randn(shape_sr, device=device)
    L = 1.0
    for _ in range(iters):
        w = adjoint_At(forward_A(v, k, sr), k, sr)
        L = float((w.norm() / (v.norm() + 1e-12)).item())
        v = w / (w.norm() + 1e-12)
    return L


def fista(y, k, sr, lam=0.01, n_iter=150, L=None, nonneg=True):
    device = y.device
    H, W = y.shape
    srH, srW = H * sr, W * sr
    if L is None:
        L = lipschitz(k, sr, (srH, srW), device)
    step = 1.0 / L
    thr = lam * step
    x = torch.zeros((srH, srW), device=device)
    z = x.clone()
    t = 1.0
    for _ in range(n_iter):
        grad = adjoint_At(forward_A(z, k, sr) - y, k, sr)
        zg = z - step * grad
        xn = F.relu(zg - thr) if nonneg else torch.sign(zg) * F.relu(zg.abs() - thr)
        tn = (1 + np.sqrt(1 + 4 * t * t)) / 2
        z = xn + ((t - 1) / tn) * (xn - x)
        x, t = xn, tn
    return x


def deconv_localize(img_np, k, sr, lam=0.01, n_iter=150, thr_frac=0.10,
                    min_sep_sr=None, device="cpu", L=None):
    """Run FISTA on a native envelope patch and return recovered peaks (SR idx) + amps + map."""
    y = torch.tensor(img_np, dtype=torch.float32, device=device)
    x_sr = fista(y, k, sr, lam=lam, n_iter=n_iter, L=L).cpu().numpy()
    if min_sep_sr is None:
        min_sep_sr = max(2, int(round(0.5 * sr)))   # ~half a native px separation
    thr = thr_frac * float(x_sr.max()) if x_sr.max() > 0 else np.inf
    peaks, amps = extract_peaks_sr(x_sr, thr, min_sep_sr=min_sep_sr)
    return peaks, amps, x_sr


def validate(densities=(3, 10, 25, 50, 80), n_each=8, lam=0.01, n_iter=150,
             thr_frac=0.10, tol_px=1.5, device="cpu", seed=3):
    sx_px, sz_px, meta = load_psf()
    k = render_psf_sr(sx_px, sz_px, 8, device=device)
    cfg0 = SynthConfig()
    sr = cfg0.sr
    # one Lipschitz for the fixed patch size (reuse across all solves)
    L = lipschitz(k, sr, (cfg0.H * sr, cfg0.W * sr), device)
    rng = np.random.default_rng(seed)
    print(f"[deconv] PSF sigma px lat {sx_px:.2f} ax {sz_px:.2f}; sr={sr}; lam={lam}; "
          f"n_iter={n_iter}; L={L:.3g}; device={device}")
    print(f"[deconv] {'dens':>5} {'GTbub':>6} {'recov':>6} {'prec':>5} {'recall':>6} "
          f"{'F1':>5} {'cnt_ratio':>9}  {'s/solve':>7}")
    rows = []
    for dens in densities:
        cfg = SynthConfig(**{**cfg0.__dict__, "dens_min": dens, "dens_max": dens})
        nm = nrec = ngt = 0
        t0 = time.time()
        for _ in range(n_each):
            s = make_sample(cfg, sx_px, sz_px, rng)
            peaks, amps, _ = deconv_localize(s["img"], k, sr, lam=lam, n_iter=n_iter,
                                             thr_frac=thr_frac, device=device, L=L)
            rec_zx = peaks / sr                          # SR idx -> native px
            gt_zx = s["locs"][:, :2]                     # (z_native, x_native)
            m, r, g = match_localizations(rec_zx, gt_zx, tol_px)
            nm += m; nrec += r; ngt += g
        dt = (time.time() - t0) / n_each
        prec = nm / nrec if nrec else 0.0
        rec = nm / ngt if ngt else 0.0
        f1 = 2 * prec * rec / (prec + rec) if (prec + rec) else 0.0
        rows.append(dict(dens=dens, gt=ngt, rec=nrec, prec=prec, recall=rec, f1=f1,
                         cnt_ratio=nrec / ngt if ngt else 0.0, s_solve=dt))
        print(f"[deconv] {dens:>5} {ngt:>6} {nrec:>6} {prec:>5.2f} {rec:>6.2f} "
              f"{f1:>5.2f} {nrec/ngt if ngt else 0:>9.2f}  {dt:>7.2f}")
    out = os.path.normpath(os.path.join(HERE, "..", "artifacts", "deconv_validation.json"))
    json.dump(rows, open(out, "w"), indent=2)
    print(f"[deconv] saved -> {out}")
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", action="store_true")
    ap.add_argument("--n_iter", type=int, default=150)
    ap.add_argument("--lam", type=float, default=0.01)
    ap.add_argument("--thr_frac", type=float, default=0.10)
    ap.add_argument("--n_each", type=int, default=8)
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    args = ap.parse_args()
    if args.validate:
        validate(n_each=args.n_each, lam=args.lam, n_iter=args.n_iter,
                 thr_frac=args.thr_frac, device=args.device)
