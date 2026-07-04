#!/usr/bin/env python
# SPDX-License-Identifier: MIT
"""QC diagnostic for the trained checkpoints: are LISTA / Deep-ULM usable, and is Deep-ULM
salvageable by re-tuning peak extraction (it has high recall but a too-dense map)?"""
import os, sys
import numpy as np
import torch
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from datagen.synth import load_psf, make_sample, SynthConfig
from deepulm.model import DeepULM
from lista.model import DeepUnfoldedULM
from common.peaks import extract_peaks_sr, match_localizations

ZART = r"Z:\Eli\phase3\artifacts"
dev = "cuda" if torch.cuda.is_available() else "cpu"
sx, sz, meta = load_psf(); cfg = SynthConfig(); sr = cfg.sr


def load(name, cls):
    ck = torch.load(os.path.join(ZART, f"{name}_ckpt.pt"), map_location=dev)
    net = cls().to(dev); net.load_state_dict(ck["state"]); net.eval()
    return net


def run(net, img):
    with torch.no_grad():
        return net(torch.from_numpy(img[None, None]).to(dev))[0, 0].cpu().numpy()


def pr(rec_zx, gt, tol=1.5):
    m, r, g = match_localizations(rec_zx, gt, tol)
    return (m / r if r else 0.0), (m / g if g else 0.0), r


deep = load("deepulm", DeepULM)
lista = load("lista", DeepUnfoldedULM)
rng = np.random.default_rng(1)

print(f"device {dev}")
for dens in [3, 25, 80]:
    c = SynthConfig(**{**cfg.__dict__, "dens_min": dens, "dens_max": dens})
    s = make_sample(c, sx, sz, rng)
    gt = s["locs"][:, :2]
    md = run(deep, s["img"]); ml = run(lista, s["img"])
    print(f"\n--- density {dens} (GT {s['n']} bubbles) ---")
    print(f"  Deep-ULM map: max {md.max():.3e} mean {md.mean():.3e} frac>0.1*max {np.mean(md>0.1*md.max()):.4f}")
    print(f"  LISTA    map: max {ml.max():.3e} mean {ml.mean():.3e}")
    print("  Deep-ULM peak sweep (thr_frac x min_sep_sr -> prec/recall/cnt):")
    for tf in [0.1, 0.3, 0.5, 0.7]:
        row = []
        for ms in [4, 8, 16]:
            thr = tf * md.max()
            pk, _ = extract_peaks_sr(md, thr, min_sep_sr=ms)
            p, r, nrec = pr(pk / sr, gt)
            row.append(f"ms{ms}:{p:.2f}/{r:.2f}/{nrec/max(1,len(gt)):.1f}")
        print(f"    thr {tf}:  " + "   ".join(row))

# visualize density-25 maps
try:
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    c = SynthConfig(**{**cfg.__dict__, "dens_min": 25, "dens_max": 25})
    s = make_sample(c, sx, sz, rng)
    md = run(deep, s["img"]); ml = run(lista, s["img"])
    fig, ax = plt.subplots(1, 4, figsize=(16, 4))
    ax[0].imshow(s["img"], cmap="magma", aspect="auto"); ax[0].set_title("input env")
    ax[1].imshow(s["target"], cmap="magma", aspect="auto"); ax[1].set_title("GT target")
    ax[2].imshow(md, cmap="magma", aspect="auto"); ax[2].set_title(f"Deep-ULM (max {md.max():.2f})")
    ax[3].imshow(ml, cmap="magma", aspect="auto"); ax[3].set_title(f"LISTA (max {ml.max():.1e})")
    fig.tight_layout(); fp = os.path.join(HERE, "artifacts", "diagnose_nets.png")
    fig.savefig(fp, dpi=110); print(f"\nsaved viz -> {fp}")
except Exception as e:
    print(f"viz skipped: {e}")
