#!/usr/bin/env python
# SPDX-License-Identifier: MIT
"""
Shared trainer for the two learned localizers (Deep-ULM, LISTA) on the ONE synthetic dataset.

Deep-STORM loss (van Sloun 2021 eq 2):  c = || f(x) - target ||^2 + zeta * || f(x) ||_1 ,  zeta=0.01
  target = Gaussian-smoothed (sigma 1 SR-px) GT delta map (the generator already produces this).
Adam lr 1e-3, on-line synthesis. Both nets train on the SAME dataset (only the model differs).

Collapse watch (Yonatan caveat): start with plain MSE+L1. If the net collapses to ~0 output, escalate
to a foreground-weighted mask loss (--mask_weight > 1; foreground = target>eps). OFF by default.

  >>> Trains on the STOPGAP PSF (phase3/artifacts/psf_training.npz). See UPGRADE_LATER.md. <<<

Usage:
  python train.py --model deepulm --iters 20000 --batch 64
  python train.py --model lista   --iters 20000 --batch 64
"""
import os, sys, time, argparse, json
import numpy as np
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from datagen.synth import load_psf, make_sample, SynthConfig                 # noqa: E402
from common.peaks import extract_peaks_sr, match_localizations              # noqa: E402
from deepulm.model import DeepULM                                            # noqa: E402
from lista.model import DeepUnfoldedULM                                      # noqa: E402

ART = os.path.join(HERE, "artifacts")


class SynthDataset(Dataset):
    """Top-level (picklable) online-synthesis dataset for DataLoader workers."""
    def __init__(self, cfg, sx_px, sz_px, length, seed=0):
        self.cfg, self.sx, self.sz, self.length, self.seed = cfg, sx_px, sz_px, length, seed

    def __len__(self):
        return self.length

    def __getitem__(self, idx):
        rng = np.random.default_rng(self.seed * 1_000_003 + idx)
        s = make_sample(self.cfg, self.sx, self.sz, rng)
        return (torch.from_numpy(s["img"][None]), torch.from_numpy(s["target"][None]))


def build_model(name):
    return DeepULM() if name == "deepulm" else DeepUnfoldedULM()


def deepstorm_loss(pred, target, zeta=0.01, mask_weight=1.0):
    if mask_weight > 1.0:
        w = 1.0 + (mask_weight - 1.0) * (target > 1e-4).float()   # weight foreground higher
        mse = (w * (pred - target) ** 2).mean()
    else:
        mse = ((pred - target) ** 2).mean()
    return mse + zeta * pred.abs().mean(), mse


@torch.no_grad()
def evaluate(net, cfg, sx, sz, device, densities=(3, 10, 25, 50, 80), n_each=6,
             thr_frac=0.10, tol_px=1.5, seed=999):
    net.eval()
    rng = np.random.default_rng(seed)
    rows = []
    for dens in densities:
        c = SynthConfig(**{**cfg.__dict__, "dens_min": dens, "dens_max": dens})
        nm = nrec = ngt = 0
        for _ in range(n_each):
            s = make_sample(c, sx, sz, rng)
            x = torch.from_numpy(s["img"][None, None]).to(device)
            m = net(x)[0, 0].cpu().numpy()
            thr = thr_frac * float(m.max()) if m.max() > 0 else np.inf
            peaks, _ = extract_peaks_sr(m, thr, min_sep_sr=max(2, int(0.5 * cfg.sr)))
            mm, r, g = match_localizations(peaks / cfg.sr, s["locs"][:, :2], tol_px)
            nm += mm; nrec += r; ngt += g
        prec = nm / nrec if nrec else 0.0
        rec = nm / ngt if ngt else 0.0
        rows.append(dict(dens=dens, gt=ngt, rec=nrec, prec=round(prec, 2), recall=round(rec, 2),
                         cnt_ratio=round(nrec / ngt, 2) if ngt else 0.0))
    net.train()
    return rows


def train(model_name, iters=20000, batch=64, lr=1e-3, zeta=0.01, mask_weight=1.0,
          device="cuda", workers=4, log_every=200, seed=0, save_every=1000, resume=False):
    cfg = SynthConfig()
    sx, sz, meta = load_psf()
    net = build_model(model_name).to(device)
    nparm = sum(p.numel() for p in net.parameters())
    opt = torch.optim.Adam(net.parameters(), lr=lr)
    ckpt_path = os.path.join(ART, f"{model_name}_ckpt.pt")
    start_it = 0
    if resume and os.path.exists(ckpt_path):
        ck = torch.load(ckpt_path, map_location=device)
        net.load_state_dict(ck["state"])
        if ck.get("opt") is not None:
            opt.load_state_dict(ck["opt"])
        start_it = int(ck.get("it", 0))
        print(f"[train] RESUMED {ckpt_path} at it {start_it}/{iters}")
    ds = SynthDataset(cfg, sx, sz, length=max(1, iters - start_it) * batch, seed=seed + start_it)
    dl = DataLoader(ds, batch_size=batch, num_workers=workers, pin_memory=True,
                    persistent_workers=(workers > 0))
    print(f"[train] {model_name}: {nparm:,} params | batch {batch} | iters {iters} | lr {lr} | "
          f"zeta {zeta} | mask_w {mask_weight} | dev {device} | save_every {save_every} | "
          f"PSF FWHM_lat {meta['final']['fwhm_lat_mm']:.3f}mm", flush=True)
    t0 = time.time()
    run = dict(model=model_name, params=int(nparm), iters=iters, batch=batch, lr=lr, zeta=zeta,
               mask_weight=mask_weight, log=[], psf=meta["final"])
    it = start_it

    def save_ckpt(at_it):
        torch.save(dict(model=model_name, state=net.state_dict(), opt=opt.state_dict(),
                        it=at_it, cfg=cfg.__dict__, psf=meta["final"]), ckpt_path)
    for x, y in dl:
        x, y = x.to(device, non_blocking=True), y.to(device, non_blocking=True)
        opt.zero_grad()
        pred = net(x)
        loss, mse = deepstorm_loss(pred, y, zeta, mask_weight)
        loss.backward(); opt.step()
        it += 1
        if it % log_every == 0 or it == 1:
            with torch.no_grad():
                pm, px = float(pred.mean()), float(pred.max())
                ym, yx = float(y.mean()), float(y.max())
            collapse = px < 0.2 * yx          # output peak far below target peak -> collapsing
            rate = it * batch / (time.time() - t0)
            lv, mv = float(loss.detach()), float(mse.detach())
            msg = (f"[train] it {it:>6}/{iters} loss {lv:.4e} mse {mv:.4e} "
                   f"pred(mean/max) {pm:.3e}/{px:.3e} tgt {ym:.3e}/{yx:.3e} "
                   f"{rate:.0f} smp/s{'  <<COLLAPSE?' if collapse else ''}")
            print(msg, flush=True)
            run["log"].append(dict(it=it, loss=lv, mse=mv, pred_max=px,
                                   tgt_max=yx, collapse=bool(collapse)))
        if it % save_every == 0:
            save_ckpt(it)                       # periodic checkpoint (resume-safe)
        if it >= iters:
            break
    # eval + save
    ev = evaluate(net, cfg, sx, sz, device)
    run["eval"] = ev
    print(f"[train] eval (recall/precision/count-ratio vs density):")
    for r in ev:
        print(f"   dens {r['dens']:>3}: prec {r['prec']:.2f} recall {r['recall']:.2f} cnt_ratio {r['cnt_ratio']:.2f}")
    save_ckpt(iters)
    json.dump(run, open(os.path.join(ART, f"{model_name}_trainlog.json"), "w"), indent=2)
    print(f"[train] saved -> {ckpt_path}  ({time.time()-t0:.0f}s total)")
    return run


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", choices=["deepulm", "lista"], required=True)
    ap.add_argument("--iters", type=int, default=20000)
    ap.add_argument("--batch", type=int, default=64)
    ap.add_argument("--lr", type=float, default=1e-3)
    ap.add_argument("--zeta", type=float, default=0.01)
    ap.add_argument("--mask_weight", type=float, default=1.0)
    ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--save_every", type=int, default=1000)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    args = ap.parse_args()
    train(args.model, iters=args.iters, batch=args.batch, lr=args.lr, zeta=args.zeta,
          mask_weight=args.mask_weight, device=args.device, workers=args.workers,
          save_every=args.save_every, resume=args.resume)
