#!/usr/bin/env python
# SPDX-License-Identifier: MIT
"""
Phase 3 sweep (NEW SPEC, config-driven via datasets.py; RESUMABLE). Runs 3 stages for one dataset:

  Stage 1  MATLAB  export_svd_envelopes(POL cache, domain)  -> env/<tag>_<domain>_env.mat
  Stage 2  Python  infer_localizers (method x tag x domain) -> locs/<method>_<tag>_<domain>_locs.mat
  Stage 3  Python  count_locrate                            -> out/locrate_<ds>.npy

NO tracking, NO QC. Two DISTINCT thresholds (do not conflate):
  clip_thr = per-domain thrFixed (the ENVELOPE clip-high threshold removing the clutter floor
             before the net; RESPONSE.md Sections 2, 5, 6). Read from readout_<ds>.mat when
             present, else passed via --thr.
  out_thr  = the SR-map PEAK threshold (which output-map peaks count as a detection). A single
             FIXED ABSOLUTE value per (method, domain), reused across all that method/domain's
             blocks. By default it is Bg-CALIBRATED: the FIRST bgFlow block is inferred with
             --bg_env (~--target_bg_rate locs/frame on Bg) and the chosen out_thr is read back
             from its locs .mat meta.out_thr, then applied (--out_thr <value>) to every remaining
             block. Supplying --out_thr overrides the calibration. The per-frame RELATIVE
             threshold is NEVER used here (infer_localizers.py RAISES without --out_thr/--bg_env).

Skips any output that already exists (resumable). Variable blocks-per-rung comes from the config.

  python run_sweep.py --dataset apr17 --domains PI --thr 409.29 --smoke
  python run_sweep.py --dataset apr17 --domains PI,fundamental,singlepol --thr 409.29,..,..
"""
import os, sys, subprocess, argparse, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import datasets as DS                                            # noqa: E402


def run_matlab(expr, matlab_exe):
    subprocess.run([matlab_exe, "-batch", expr], check=True)


def read_out_thr(locs_file):
    """Read the calibrated FIXED ABSOLUTE out_thr back from an infer output .mat (meta.out_thr).
    Raises if the file was not produced by a --bg_env / --out_thr run (meta.out_thr == -1 sentinel)."""
    from scipy.io import loadmat
    meta = loadmat(locs_file, squeeze_me=True, struct_as_record=False)["meta"]
    val = float(np.atleast_1d(meta.out_thr).ravel()[0])
    if not np.isfinite(val) or val < 0:
        raise SystemExit("[sweep] %s has no calibrated out_thr (meta.out_thr=%r); it was not "
                         "produced by a --bg_env calibration (or --out_thr) run." % (locs_file, val))
    return val


def read_thrFixed(readout_file, domains):
    """Per-domain thrFixed from readout_<ds>.mat (var readout). Returns {domain: thr} or None.
    Reuses count_locrate.read_anchor_readout, which handles BOTH -v7 (scipy) and -v7.3/HDF5
    (the real RUNNER writes readout as -v7.3, which scipy.loadmat cannot read)."""
    if not (readout_file and os.path.exists(readout_file)):
        return None
    import count_locrate as CL
    ro = CL.read_anchor_readout(readout_file)
    rd_dom = list(ro["domains"])
    thr = np.asarray(ro["thrFixed"], float).ravel()
    return {dom: float(thr[rd_dom.index(dom)]) for dom in domains if dom in rd_dom}


def resolve_thresholds(cfg, domains, thr_arg):
    """thrFixed per domain: readout_<ds>.mat if present, else --thr (comma list matching --domains)."""
    thr = read_thrFixed(cfg.get("readoutFile"), domains)
    if thr and all(d in thr for d in domains):
        print("[sweep] thrFixed from %s: %s" % (cfg["readoutFile"], thr))
        return thr
    if not thr_arg:
        raise SystemExit("[sweep] readout_%s.mat absent and no --thr given; pass per-domain "
                         "thrFixed as --thr (comma list matching --domains)." % cfg["name"])
    vals = [float(v) for v in str(thr_arg).split(",")]
    if len(vals) == 1:
        vals = vals * len(domains)
    if len(vals) != len(domains):
        raise SystemExit("[sweep] --thr count (%d) must match --domains count (%d)" % (len(vals), len(domains)))
    out = dict(zip(domains, vals))
    print("[sweep] thrFixed from --thr: %s" % out)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", required=True, choices=list(DS.DATASETS))
    ap.add_argument("--domains", default="PI", help="comma list (default PI; allow all 3)")
    ap.add_argument("--methods", default="deepulm,lista,deconv")
    ap.add_argument("--root", default=None, help="scratch root (default <outDir>\\infer_run)")
    ap.add_argument("--thr", default=None, help="per-domain thrFixed (comma list matching --domains)")
    ap.add_argument("--out_thr", default=None, help="per-domain fixed ABSOLUTE SR-map threshold (comma list); "
                    "overrides the per-(method,domain) Bg calibration")
    ap.add_argument("--target_bg_rate", type=float, default=20.0,
                    help="target Bg loc-rate/frame for the per-(method,domain) out_thr calibration (default 20)")
    ap.add_argument("--ckptdir", default=r"Q:\Eli\phase3\artifacts")
    ap.add_argument("--matlab", default=r"C:\Program Files\MATLAB\R2025b\bin\matlab.exe")
    ap.add_argument("--frames", type=int, default=0, help="export-frame cap (0=all); --smoke sets 400")
    ap.add_argument("--deconv_frames", type=int, default=400, help="frame subset for the slow FISTA deconv")
    ap.add_argument("--use_gpu", type=int, default=1)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--no_bg_static", action="store_true", help="skip bgStatic blocks")
    args = ap.parse_args()

    cfg = DS.get_dataset(args.dataset)
    domains = args.domains.split(",")
    methods = args.methods.split(",")
    thr = resolve_thresholds(cfg, domains, args.thr)
    out_thr = None
    if args.out_thr:
        vals = [float(v) for v in args.out_thr.split(",")]
        if len(vals) == 1:
            vals = vals * len(domains)
        out_thr = dict(zip(domains, vals))

    root = args.root or os.path.join(cfg["outDir"], "infer_run")
    envdir = os.path.join(root, "env"); locdir = os.path.join(root, "locs"); outdir = os.path.join(root, "out")
    for d in (envdir, locdir, outdir):
        os.makedirs(d, exist_ok=True)

    tags = DS.all_block_tags(cfg, include_bg_static=not args.no_bg_static)
    frames = args.frames
    if args.smoke:
        rb = DS.rung_block_tags(cfg)
        tags = [rb[0], rb[len(rb) // 2], cfg["bgFlow"][0]]
        frames = frames or 400
    nFr_ml = "[]" if frames == 0 else str(frames)
    gpu_ml = "true" if args.use_gpu else "false"
    print("[sweep] %s | domains=%s | methods=%s | tags=%d | smoke=%s | export frames=%s | deconv frames=%d" % (
        cfg["name"], domains, methods, len(tags), args.smoke, nFr_ml, args.deconv_frames))

    # ---- Stage 1: export envelopes per (tag, domain) from the POL caches (skip existing) ----
    todo = []
    for t in tags:
        for dom in domains:
            envp = os.path.join(envdir, "%s_%s_env.mat" % (t, dom))
            if not os.path.exists(envp):
                todo.append((t, dom, envp))
    if todo:
        print("[sweep] Stage 1: export %d envelopes (MATLAB useGPU=%s)" % (len(todo), gpu_ml))
        # MATLAB single-quoted strings treat backslash literally -> pass Windows paths AS-IS
        # (do NOT escape backslashes). The whole expr is one argv to matlab.exe (no shell).
        polcell = "{" + ",".join("'%s'" % DS.pol_cache_path(cfg, t) for t, _, _ in todo) + "}"
        domcell = "{" + ",".join("'%s'" % d for _, d, _ in todo) + "}"
        outcell = "{" + ",".join("'%s'" % p for _, _, p in todo) + "}"
        expr = ("cd('%s'); pol=%s; dom=%s; outp=%s; "
                "for i=1:numel(pol), assert(isfile(pol{i}),'no POL cache %%s',pol{i}); "
                "export_svd_envelopes(pol{i}, outp{i}, dom{i}, %s, %s); end" % (
                    HERE, polcell, domcell, outcell, nFr_ml, gpu_ml))
        run_matlab(expr, args.matlab)
    else:
        print("[sweep] Stage 1: all envelopes present, skipping")

    INFER = os.path.join(HERE, "infer_localizers.py")

    # ---- Stage 2: infer per (method, tag, domain) (skip existing, resumable). clip_thr=thrFixed
    #      per domain (unchanged). out_thr is ONE FIXED ABSOLUTE value per (method, domain): either
    #      supplied via --out_thr, or Bg-CALIBRATED by inferring the FIRST bgFlow block with
    #      --bg_env (~--target_bg_rate/frame) and reading meta.out_thr back, then reused for every
    #      remaining block. The per-frame RELATIVE threshold is NEVER passed here. ----
    for m in methods:
        ck = ["--ckpt", os.path.join(args.ckptdir, "%s_ckpt.pt" % m)] if m in ("deepulm", "lista") else []
        nf = ["--n_frames", str(args.deconv_frames)] if m == "deconv" else []
        for dom in domains:
            cthr = ["--clip_thr", str(thr[dom])]
            bg0 = cfg["bgFlow"][0]
            bg0_env = os.path.join(envdir, "%s_%s_env.mat" % (bg0, dom))
            bg0_out = os.path.join(locdir, "%s_%s_%s_locs.mat" % (m, bg0, dom))

            if out_thr and dom in out_thr:
                # explicit per-domain absolute override -> no calibration; bg0 is inferred in-loop
                othr_val = out_thr[dom]
                calibrated_bg0 = False
                print("[sweep] %s/%s: out_thr=%.6g (supplied --out_thr)" % (m, dom, othr_val))
            else:
                # Bg-calibrate on the FIRST bgFlow block (produces its locs + meta.out_thr)
                calibrated_bg0 = True
                if not os.path.exists(bg0_out):
                    if not os.path.exists(bg0_env):
                        raise SystemExit("[sweep] %s/%s: cannot calibrate out_thr -- missing Bg "
                                         "envelope %s (run Stage 1 first)." % (m, dom, bg0_env))
                    t0 = time.time()
                    print("[sweep] infer %s %s %s (Bg-calibrate out_thr @ %.3g/frame) ..." % (
                        m, bg0, dom, args.target_bg_rate), flush=True)
                    subprocess.run([sys.executable, INFER, "--method", m, "--env", bg0_env,
                                    "--out", bg0_out, "--domain", dom, "--bg_env", bg0_env,
                                    "--target_bg_rate", str(args.target_bg_rate)] + cthr + ck + nf,
                                   check=True)
                    print("[sweep]   (%.0fs)" % (time.time() - t0))
                othr_val = read_out_thr(bg0_out)
                print("[sweep] %s/%s: out_thr=%.6g (Bg-calibrated from %s)" % (
                    m, dom, othr_val, os.path.basename(bg0_out)))

            othr = ["--out_thr", str(othr_val)]
            for t in tags:
                if calibrated_bg0 and t == bg0:
                    continue                                  # already produced during calibration
                envp = os.path.join(envdir, "%s_%s_env.mat" % (t, dom))
                out = os.path.join(locdir, "%s_%s_%s_locs.mat" % (m, t, dom))
                if os.path.exists(out):
                    continue
                if not os.path.exists(envp):
                    print("[sweep]   skip %s/%s/%s: no env" % (m, t, dom)); continue
                t0 = time.time(); print("[sweep] infer %s %s %s ..." % (m, t, dom), flush=True)
                subprocess.run([sys.executable, INFER, "--method", m, "--env", envp, "--out", out,
                                "--domain", dom] + cthr + othr + ck + nf, check=True)
                print("[sweep]   (%.0fs)" % (time.time() - t0))

    # ---- Stage 3: count (Python; needs the full rung set, skip in smoke) ----
    if args.smoke:
        print("[sweep] SMOKE OK: export + infer ran. (Stage 3 count needs the full rung set -> "
              "run without --smoke.)")
        return
    print("[sweep] Stage 3: count_locrate")
    import count_locrate
    outp = os.path.join(outdir, "locrate_%s.npy" % cfg["name"])
    count_locrate.run_dataset(cfg, locdir, out_path=outp, domains=domains, methods=tuple(methods),
                              roi_file=cfg["roiFile"])
    print("[sweep] DONE -> %s" % outp)


if __name__ == "__main__":
    main()
