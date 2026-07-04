#!/usr/bin/env python
# SPDX-License-Identifier: MIT
"""
Phase 3 inference (NEW SPEC, REANALYSIS_INPUTS_RESPONSE.md Sections 2-4): real SVD-filtered
envelope -> NATIVE-GRID PIXEL localizations [row, col, amp, frameIdx] for a chosen localizer
(deconv | deepulm | lista). The pixel [row, col] index is what the new in-tube-localization-rate
count (count_locrate.py) needs for its ROI-mask lookup; there is NO tracking and NO QC anymore.

Per-block grid: the grid `g` is loaded from the exported envelope .mat itself (each block carries
its own natural-depth grid, Section 4a) -- NOT from a global batch_config.

CLIP-HIGH preprocessing (validated 2026-06-26): the synthetic-trained nets FLOOD on the real SVD
envelope because it has a high residual-clutter floor they never saw. Clipping the input at the
per-domain fixed detection threshold (== readout.thrFixed for that domain, passed via --clip_thr)
removes that floor so the real input matches the clean-background synthetic training domain. All
learned methods share this single per-domain threshold (= the classical detector's operating point)
-> fair comparison. The nets / clip-high run on abs(IQf) (Section 2 step 3) -- that is what `env`
already is (export_svd_envelopes.m saved env = single(abs(IQf))).

Output peak threshold is a FIXED ABSOLUTE SR-map threshold, reused across all blocks of a
(method, domain): either supplied via --out_thr, or Bg-calibrated via --bg_env <bg_env.mat>
(calibrate_out_thr runs the net on the Bg envelope and returns the absolute threshold giving
~--target_bg_rate localizations/frame on Bg; saved to meta.out_thr so the sweep can read it
back). A fixed absolute threshold keeps the count monotonic in concentration. The per-frame
RELATIVE threshold (out_thr_frac*map.max) FLOODS bubble-free frames and is a DEBUG-only opt-in
via explicit --out_thr_frac; with none of --out_thr / --bg_env / --out_thr_frac given this script
RAISES (it never silently uses the relative path). The connected-component peak extractor
(common/peaks.py) is O(N) and cannot hang.
"""
import os, sys, argparse, time
import numpy as np
import h5py
import torch
from scipy.io import savemat

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
from common.peaks import extract_peaks_sr                          # noqa: E402
from deepulm.model import DeepULM                                  # noqa: E402
from lista.model import DeepUnfoldedULM                            # noqa: E402
from deconv.fista_deconv import render_psf_sr, fista, lipschitz    # noqa: E402
from datagen.synth import load_psf                                 # noqa: E402

CLIP_SCALE = 200.0            # maps above-threshold bubble signal to ~O(1) (synthetic training scale)
SR = 8


def load_env(path):
    with h5py.File(path, "r") as f:
        e = f["env"][()]                                  # (nFr, nX, nZ) HDF5-reversed of [nZ,nX,nFr]
    return np.transpose(e, (0, 2, 1)).astype(np.float32)  # [nFr, nZ, nX]


def _rd(group, name):
    return np.asarray(group[name][()]).ravel()


def load_grid_from_env(path):
    """Per-block grid g from the exported envelope .mat (Section 4a): g.xGrid/zGrid/dx/dz/nX/nZ
    (+ minSep_mm). MATLAB -v7.3 stores the struct as an HDF5 group; read its datasets directly.
    Returns (grid_dict, nFrames)."""
    with h5py.File(path, "r") as f:
        g = f["g"]
        grid = dict(
            xGrid=_rd(g, "xGrid").astype(np.float64),
            zGrid=_rd(g, "zGrid").astype(np.float64),
            dx=float(_rd(g, "dx")[0]),
            dz=float(_rd(g, "dz")[0]),
            nX=int(_rd(g, "nX")[0]),
            nZ=int(_rd(g, "nZ")[0]),
        )
        grid["minSep_mm"] = float(_rd(g, "minSep_mm")[0]) if "minSep_mm" in g else float("nan")
        nFrames = int(np.asarray(f["nFrames"][()]).ravel()[0]) if "nFrames" in f else None
    return grid, nFrames


def preprocess(env, clip_thr, clip_scale=CLIP_SCALE):
    """Clip-high at the per-domain fixed threshold then rescale (Section 2 step 3 domain)."""
    return (np.maximum(env - clip_thr, 0.0) / clip_scale).astype(np.float32)


def sr_to_pixel(peaks_sr, sr, nZ, nX):
    """SR-grid peak (z_sr, x_sr) -> NATIVE-grid pixel index [row, col] on g (Section 2 step 3c).
    row = round(z_sr/sr) clipped to [0,nZ-1]; col = round(x_sr/sr) clipped to [0,nX-1]."""
    row = np.clip(np.round(peaks_sr[:, 0] / sr).astype(np.int64), 0, nZ - 1)
    col = np.clip(np.round(peaks_sr[:, 1] / sr).astype(np.int64), 0, nX - 1)
    return row, col


# ---------- per-frame map producers ----------
def make_map_fn(method, device, ckpt=None, lam=0.10, n_iter=120):
    """Return a function frames[B,H,W] -> list of B SR maps, and the SR factor.
    Deconv lambda defaults to the FISTA spec value 0.10 (RESPONSE.md Section 5)."""
    if method in ("deepulm", "lista"):
        cls = DeepULM if method == "deepulm" else DeepUnfoldedULM
        ck = torch.load(ckpt, map_location=device)
        net = cls().to(device); net.load_state_dict(ck["state"]); net.eval()

        def fn(frames):
            with torch.no_grad():
                x = torch.from_numpy(frames[:, None]).to(device)
                return [m for m in net(x)[:, 0].cpu().numpy()]
        return fn, SR

    elif method == "deconv":
        sx, sz, _ = load_psf()
        k = render_psf_sr(sx, sz, SR, device=device)
        state = {"L": None}

        def fn(frames):
            H, W = frames.shape[1], frames.shape[2]
            if state["L"] is None:
                state["L"] = lipschitz(k, SR, (H * SR, W * SR), device)
            out = []
            for f in frames:
                y = torch.from_numpy(f).to(device)
                out.append(fista(y, k, SR, lam=lam, n_iter=n_iter, L=state["L"]).cpu().numpy())
            return out
        return fn, SR
    raise ValueError(method)


def maps_to_locs(map_fn, sr, envc, nZ, nX, out_thr_abs, batch, min_sep_sr,
                 out_thr_frac=0.2, frame_offset=0):
    """Per frame: net/deconv SR map -> peaks -> NATIVE pixel [row, col] -> [row, col, amp, frame].
    out_thr_abs (fixed absolute SR-map threshold) is the operating mode for the real sweep; the
    relative fallback (out_thr_frac*map.max) is opt-in only and never used by run_sweep.py."""
    locs = []
    for i in range(0, len(envc), batch):
        for j, M in enumerate(map_fn(envc[i:i + batch])):
            fr = i + j
            mx = float(M.max())
            if mx <= 0:
                continue
            thr = out_thr_abs if out_thr_abs is not None else out_thr_frac * mx
            pk, amp = extract_peaks_sr(M, thr, min_sep_sr)
            if len(pk):
                row, col = sr_to_pixel(pk, sr, nZ, nX)
                locs.append(np.column_stack(
                    [row, col, amp, np.full(len(pk), fr + 1 + frame_offset, np.float32)]))
    return np.vstack(locs).astype(np.float32) if locs else np.zeros((0, 4), np.float32)


def calibrate_out_thr(map_fn, sr, bg_envc, target_rate, min_sep_sr, batch):
    """Bg-calibrate the FIXED ABSOLUTE SR-map peak threshold (FIX 1).

    Runs `map_fn` over the (already clip-high preprocessed) Bg envelope `bg_envc` and returns the
    single absolute out_thr whose strict-greater peak count yields ~`target_rate` localizations per
    frame on Bg. With no tracking to reject clutter, this fixed Bg-referenced operating point is
    what keeps the count monotonic in concentration (a per-frame relative threshold floods
    bubble-free Bg). It is NON-circular: tuned to a Bg false-alarm rate, never to the anchor count.

    Candidate peaks are collected above a SMALL per-map floor `1e-3 * map.max` (NOT ~0: a near-zero
    floor lights up every flat pixel and floods) via the O(N) connected-component extractor
    `extract_peaks_sr` (no O(N^2) loop -> cannot hang). Because `extract_peaks_sr`'s local-maximum
    set does NOT depend on the threshold (only the `> thr` filter does), the candidates above the
    floor are a superset of every peak any `thr >= floor` would return -- so the K-th largest
    candidate amplitude (K = round(target_rate * nFr)) is the EXACT absolute threshold that yields
    ~K Bg peaks total (==> ~target_rate / frame). `sr` is accepted for signature parity with the
    map producers and is not needed here.
    """
    nfr = int(len(bg_envc))
    if nfr == 0:
        raise ValueError("calibrate_out_thr: empty Bg envelope")
    amps_all = []
    for i in range(0, nfr, batch):
        for M in map_fn(bg_envc[i:i + batch]):
            mx = float(M.max())
            if mx <= 0:
                continue
            floor = 1e-3 * mx                                  # small per-map floor (NOT ~0)
            _pk, amp = extract_peaks_sr(M, floor, min_sep_sr)  # O(N) connected components
            if len(amp):
                amps_all.append(amp)
    if not amps_all:
        # The map producer fired on NO Bg frame (e.g. FISTA deconvolution: its L1 penalty drives the
        # clip-high'd clutter-only Bg to exactly zero). That is the IDEAL clean-background case: there
        # is no false floor to calibrate a threshold against, so the operating point is 0 -- every
        # (sparse) nonzero detection in the data blocks counts, and Bg contributes ~0 to the pedestal.
        # (extract_peaks_sr uses a strict `> thr`, so thr=0 still excludes the exact zeros.)
        print("[calibrate] no Bg peaks (self-thresholded clean Bg) -> out_thr = 0", flush=True)
        return 0.0
    amps = np.concatenate(amps_all)                            # all candidate amplitudes (>0)
    K = max(1, int(round(float(target_rate) * nfr)))           # target total Bg peaks
    amps_desc = np.sort(amps)[::-1]
    if K >= amps.size:
        # Bg cannot reach target_rate even admitting every candidate -> lowest operating point that
        # still admits them all (best effort; resulting Bg rate < target_rate).
        return float(amps_desc[-1]) * (1.0 - 1e-6)
    return float(amps_desc[K])                                 # (K+1)-th largest -> ~K strictly above


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--method", choices=["deconv", "deepulm", "lista"], required=True)
    ap.add_argument("--env", required=True, help="block SVD-envelope .mat (vars env, g, nFrames)")
    ap.add_argument("--out", required=True, help="output localizations .mat (var localizations [N x 4])")
    ap.add_argument("--domain", default="PI", choices=["PI", "fundamental", "singlepol"],
                    help="domain bookkeeping (the env was already built for this domain)")
    ap.add_argument("--ckpt", default=None, help="net checkpoint (deepulm/lista)")
    ap.add_argument("--clip_thr", type=float, required=True,
                    help="per-domain FIXED clip-high threshold == readout.thrFixed(domain)")
    ap.add_argument("--out_thr", type=float, default=None,
                    help="fixed ABSOLUTE SR-map peak threshold (operating mode for the sweep)")
    ap.add_argument("--bg_env", default=None,
                    help="Bg-block envelope .mat to CALIBRATE a fixed absolute out_thr "
                         "(~--target_bg_rate locs/frame on Bg); saved to meta.out_thr")
    ap.add_argument("--target_bg_rate", type=float, default=20.0,
                    help="target Bg loc-rate/frame for --bg_env calibration (default 20)")
    ap.add_argument("--out_thr_frac", type=float, default=None,
                    help="DEBUG ONLY relative SR-map threshold = frac*map.max; explicit opt-in, "
                         "NEVER used by the real sweep (it floods bubble-free frames)")
    ap.add_argument("--lam", type=float, default=0.10, help="FISTA lambda (deconv); spec = 0.10")
    ap.add_argument("--min_sep_sr", type=int, default=4)
    ap.add_argument("--batch", type=int, default=2)
    ap.add_argument("--n_frames", type=int, default=None)
    ap.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu")
    args = ap.parse_args()

    t0 = time.time()
    grid, nFrames_cache = load_grid_from_env(args.env)        # per-block g (Section 4a)
    nZ, nX = grid["nZ"], grid["nX"]
    env = load_env(args.env)
    if args.n_frames:
        env = env[:args.n_frames]
    nFrames = int(len(env))
    envc = preprocess(env, args.clip_thr)                     # clip-high at thrFixed (Section 2/3)
    map_fn, sr = make_map_fn(args.method, args.device, ckpt=args.ckpt, lam=args.lam)

    if args.out_thr is not None:                              # 1) supplied fixed absolute
        out_thr_abs = args.out_thr
        thr_desc = "abs %.4g (fixed, supplied)" % out_thr_abs
    elif args.bg_env is not None:                             # 2) Bg-calibrate a fixed absolute
        bg_env = load_env(args.bg_env)
        if args.n_frames:
            bg_env = bg_env[:args.n_frames]
        bg_envc = preprocess(bg_env, args.clip_thr)           # clip Bg at the SAME thrFixed
        out_thr_abs = calibrate_out_thr(map_fn, sr, bg_envc, args.target_bg_rate,
                                        args.min_sep_sr, args.batch)
        thr_desc = "abs %.4g (Bg-calibrated @ %.3g/frame on %s)" % (
            out_thr_abs, args.target_bg_rate, os.path.basename(args.bg_env))
    elif args.out_thr_frac is not None:                       # 3) DEBUG-only relative opt-in
        out_thr_abs = None
        thr_desc = "%.3g*max/frame (RELATIVE DEBUG fallback; NOT for the sweep)" % args.out_thr_frac
    else:                                                     # 4) sweep path with no threshold -> RAISE
        raise SystemExit(
            "[infer] no SR-map peak threshold given: pass --out_thr <abs> or --bg_env <bg_env.mat> "
            "(Bg-calibration). The per-frame relative threshold floods bubble-free frames and is "
            "DEBUG-only via an explicit --out_thr_frac; it is never used by the real sweep.")

    locs = maps_to_locs(map_fn, sr, envc, nZ, nX, out_thr_abs, args.batch, args.min_sep_sr,
                        out_thr_frac=(args.out_thr_frac if args.out_thr_frac is not None else 0.2))

    # grid struct saved alongside so count_locrate.py is self-contained (loads g + nFrames here)
    g_out = dict(xGrid=grid["xGrid"], zGrid=grid["zGrid"], dx=grid["dx"], dz=grid["dz"],
                 nX=nX, nZ=nZ, minSep_mm=grid["minSep_mm"])
    savemat(args.out, {
        "localizations": locs,                               # [N x 4] = [row, col, amp, frameIdx]
        "g": g_out,
        "nFrames": nFrames,
        "meta": {"method": args.method, "domain": args.domain, "clip_thr": args.clip_thr,
                 "clip_scale": CLIP_SCALE, "out_thr": (out_thr_abs if out_thr_abs is not None else -1.0),
                 "sr": sr, "nFrames": nFrames, "nZ": nZ, "nX": nX,
                 "cols": "row,col,amp,frameIdx (row=z pixel 0-based, col=x pixel 0-based)"},
    })
    rate = len(locs) / max(1, nFrames)
    print("[infer] %s/%s: %d locs over %d frames (%.1f/frame) | thr %s | clip %.4g | %.0fs -> %s" % (
        args.method, args.domain, len(locs), nFrames, rate, thr_desc, args.clip_thr,
        time.time() - t0, args.out))


if __name__ == "__main__":
    main()
