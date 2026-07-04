#!/usr/bin/env python
# SPDX-License-Identifier: MIT
"""
Phase 3 SHARED synthetic training-data generator (ONE dataset for all three localizers).

  >>> Uses the STOPGAP PSF (phase3/artifacts/psf_training.npz). See phase3/UPGRADE_LATER.md:
  when the followup measured LINEAR PSF lands, regenerate this dataset and RETRAIN both nets. <<<

Per van Sloun 2021 (Deep-ULM) Sec II.A / Sec V and the Phase-3 handoff:
  - online-synthesized point sources convolved with the measured anisotropic-Gaussian PSF
  - span the density range INCLUDING past the classical localization ceiling (the net must learn
    to resolve OVERLAPPING point-spread functions -- this is the whole point)
  - augmentation: multiplicative PSF-parameter variance ~ N(1, 0.1); additive 2% white + 5% colored
    noise (van Sloun 2021 Sec II.A: colored = white noise filtered by a 2D Gaussian, sigma 1.2 px)
  - NO real clutter in training (preserves the sim->real generalization test; real Bg blocks are
    used only at INFERENCE to calibrate the peak false-alarm rate -- the non-circular threshold)
  - emits, per sample: native-res envelope image (net/deconv INPUT), the GT bubble locations
    (sub-pixel, native coords), and the x8 super-res Deep-STORM target (Gaussian-smoothed deltas)

The same samples feed: sparse deconvolution (img + locs), Deep-ULM (img -> target), LISTA (img -> target).

Axis convention matches the data: z=axial (rows, the 404 dim), x=lateral (cols, the 193 dim).
Works in NORMALIZED units (bubble peak ~1); inference must normalize real envelope to match.
"""
import os, json
import numpy as np

HERE = os.path.dirname(__file__)
PSF_ARTIFACT = os.path.normpath(os.path.join(HERE, "..", "artifacts", "psf_training.npz"))
DX_MM = DZ_MM = 0.0346875
FWHM_K = 2.0 * np.sqrt(2.0 * np.log(2.0))


def load_psf(path=PSF_ARTIFACT):
    d = np.load(path, allow_pickle=True)
    sx_px = float(d["sigma_x_mm"]) / DX_MM      # lateral sigma in native px
    sz_px = float(d["sigma_z_mm"]) / DZ_MM      # axial sigma in native px
    meta = json.loads(str(d["meta"]))
    return sx_px, sz_px, meta


class SynthConfig:
    def __init__(self, H=64, W=64, sr=8,
                 dens_min=1, dens_max=80,            # bubbles per patch (log-uniform); max -> heavy overlap
                 amp_logsigma=0.4,                   # bubble amplitude ~ lognormal(0, amp_logsigma)
                 psf_jitter=0.10,                    # multiplicative N(1, psf_jitter) on sigma_x, sigma_z
                 white=0.02, colored=0.05,           # additive noise stds (fraction of unit amp)
                 target_sigma_sr=1.0,                # Deep-STORM target Gaussian smoothing (SR px)
                 edge_margin_px=1.0):
        self.__dict__.update(locals()); del self.__dict__["self"]


def _render_gauss(img, zc, xc, amp, sz, sx, nsig=4.0):
    """Add one anisotropic Gaussian at sub-pixel native coord (zc,xc) to img in-place."""
    H, W = img.shape
    z0, z1 = max(0, int(zc - nsig * sz)), min(H, int(zc + nsig * sz) + 1)
    x0, x1 = max(0, int(xc - nsig * sx)), min(W, int(xc + nsig * sx) + 1)
    if z0 >= z1 or x0 >= x1:
        return
    zz = np.arange(z0, z1)[:, None]; xx = np.arange(x0, x1)[None, :]
    img[z0:z1, x0:x1] += amp * np.exp(-((zz - zc) ** 2) / (2 * sz ** 2) - ((xx - xc) ** 2) / (2 * sx ** 2))


def _colored_noise(H, W, rng, sigma_px=1.2):
    """Colored background noise per van Sloun 2021 Sec II.A: white noise spatially filtered by a
    2D Gaussian (sigma 1.2 px). Returned at unit std (scaled by cfg.colored at the call site)."""
    from scipy.ndimage import gaussian_filter
    n = gaussian_filter(rng.standard_normal((H, W)).astype(np.float32), sigma=sigma_px, mode="reflect")
    s = n.std()
    return (n / s) if s > 0 else n


_TKERNEL = {}


def _target_kernel(sigma_sr):
    """Cached PEAK-normalized 2D Gaussian stamp (peak=1) so the target peak at a bubble == its
    amplitude (O(1)). Peak-normalization (not sum=1, which gave peak ~0.16) gives the net a strong
    enough positive signal to avoid the all-zero/constant Deep-STORM collapse. (fix 2026-06-25)"""
    key = round(sigma_sr, 4)
    if key not in _TKERNEL:
        r = max(1, int(np.ceil(3 * sigma_sr)))
        zz, xx = np.mgrid[-r:r + 1, -r:r + 1]
        g = np.exp(-(zz ** 2 + xx ** 2) / (2 * sigma_sr ** 2)).astype(np.float32)
        _TKERNEL[key] = (g / g.max(), r)
    return _TKERNEL[key]


def _stamp(target, zi, xi, amp, g, r):
    H, W = target.shape
    z0, z1 = max(0, zi - r), min(H, zi + r + 1)
    x0, x1 = max(0, xi - r), min(W, xi + r + 1)
    if z0 >= z1 or x0 >= x1:
        return
    gz0, gx0 = z0 - (zi - r), x0 - (xi - r)
    target[z0:z1, x0:x1] += amp * g[gz0:gz0 + (z1 - z0), gx0:gx0 + (x1 - x0)]


def make_sample(cfg, sx_px, sz_px, rng):
    """Return dict(img[H,W], target[srH,srW], locs[N,3]=(z_native,x_native,amp), n, density)."""
    H, W, sr = cfg.H, cfg.W, cfg.sr
    n = int(round(np.exp(rng.uniform(np.log(cfg.dens_min), np.log(cfg.dens_max)))))
    m = cfg.edge_margin_px
    zc = rng.uniform(m, H - 1 - m, n)
    xc = rng.uniform(m, W - 1 - m, n)
    amp = np.exp(rng.normal(0.0, cfg.amp_logsigma, n)).astype(np.float32)
    # per-patch PSF jitter (robustness to +-10% PSF / depth variation)
    sxa = sx_px * (1.0 + rng.normal(0, cfg.psf_jitter))
    sza = sz_px * (1.0 + rng.normal(0, cfg.psf_jitter))
    sxa, sza = max(sxa, 0.5), max(sza, 0.4)
    img = np.zeros((H, W), np.float32)
    for i in range(n):
        _render_gauss(img, zc[i], xc[i], amp[i], sza, sxa)
    # noise: additive white + colored (at PSF scale), then rectify (envelope is non-negative)
    img += cfg.white * rng.standard_normal((H, W)).astype(np.float32)
    img += cfg.colored * _colored_noise(H, W, rng)
    np.maximum(img, 0.0, out=img)
    # GT super-res target (Deep-STORM): Gaussian stamps on the x8 grid (fast; emulates
    # gaussian_filter(deltas) but O(n_bubbles * stamp) instead of O(srH*srW) per sample)
    target = np.zeros((H * sr, W * sr), np.float32)
    g, r = _target_kernel(cfg.target_sigma_sr)
    for i in range(n):
        _stamp(target, int(round(zc[i] * sr)), int(round(xc[i] * sr)), amp[i], g, r)
    locs = np.column_stack([zc, xc, amp]).astype(np.float32)
    return dict(img=img, target=target, locs=locs, n=n,
                density=n / (H * W), psf=(sza, sxa))


# ---- PyTorch Dataset (online synthesis) ----
def make_dataset(cfg=None, seed=0, length=20000):
    import torch
    from torch.utils.data import Dataset
    if cfg is None:
        cfg = SynthConfig()
    sx_px, sz_px, _ = load_psf()

    class _DS(Dataset):
        def __init__(self):
            self.cfg, self.sx, self.sz, self.length = cfg, sx_px, sz_px, length
        def __len__(self):
            return self.length
        def __getitem__(self, idx):
            rng = np.random.default_rng(seed * 1_000_003 + idx)   # deterministic per-index (resumable)
            s = make_sample(self.cfg, self.sx, self.sz, rng)
            x = torch.from_numpy(s["img"][None])                  # [1,H,W]
            y = torch.from_numpy(s["target"][None])               # [1,srH,srW]
            return x, y
    return _DS(), cfg


def dump_validation_set(path, n=200, cfg=None, seed=12345):
    """Fixed validation set for evaluation/QC (img, target, locs)."""
    if cfg is None:
        cfg = SynthConfig()
    sx_px, sz_px, meta = load_psf()
    rng = np.random.default_rng(seed)
    samples = [make_sample(cfg, sx_px, sz_px, rng) for _ in range(n)]
    np.savez_compressed(path,
                        imgs=np.stack([s["img"] for s in samples]),
                        targets=np.stack([s["target"] for s in samples]),
                        locs=np.array([s["locs"] for s in samples], dtype=object),
                        ns=np.array([s["n"] for s in samples]),
                        cfg=json.dumps(cfg.__dict__), psf=json.dumps(meta["final"]))
    return samples


def qc_figure(path, cfg=None):
    """Render low/mid/high-density samples + their SR targets for visual inspection."""
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    if cfg is None:
        cfg = SynthConfig()
    sx_px, sz_px, meta = load_psf()
    rng = np.random.default_rng(7)
    fig, axes = plt.subplots(2, 3, figsize=(12, 7))
    for j, dens in enumerate([3, 25, 80]):
        c = SynthConfig(**{**cfg.__dict__, "dens_min": dens, "dens_max": dens})
        s = make_sample(c, sx_px, sz_px, rng)
        axes[0, j].imshow(s["img"], cmap="magma", aspect="auto")
        axes[0, j].plot(s["locs"][:, 1], s["locs"][:, 0], "c.", ms=3)
        axes[0, j].set_title(f"input env, {s['n']} bubbles")
        axes[1, j].imshow(s["target"], cmap="magma", aspect="auto")
        axes[1, j].set_title(f"x{cfg.sr} SR target")
    fig.suptitle(f"Synthetic data (STOPGAP PSF: lat {meta['final']['fwhm_lat_mm']:.3f}mm / "
                 f"ax {meta['final']['fwhm_ax_mm']:.3f}mm) — UPGRADE_LATER.md", color="#b00", fontsize=9)
    fig.tight_layout(); fig.savefig(path, dpi=120)
    print(f"[synth] QC figure -> {path}")


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--qc", action="store_true")
    ap.add_argument("--dump", type=int, default=0)
    args = ap.parse_args()
    sx_px, sz_px, meta = load_psf()
    print(f"[synth] PSF sigma px: lateral {sx_px:.3f}  axial {sz_px:.3f}  (FWHM lat "
          f"{meta['final']['fwhm_lat_mm']:.4f} mm)")
    art = os.path.normpath(os.path.join(HERE, "..", "artifacts"))
    if args.qc:
        qc_figure(os.path.join(art, "synth_QC.png"))
    if args.dump:
        s = dump_validation_set(os.path.join(art, "synth_val.npz"), n=args.dump)
        print(f"[synth] dumped {len(s)} validation samples (mean bubbles/patch "
              f"{np.mean([x['n'] for x in s]):.1f})")
