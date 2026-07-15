# FINDINGS 2026-07-14 — detection-threshold provenance + the PSF number

Sequel to `FINDINGS_2026_06_11_lowconc_floor.md` (internal; documented the PREVIOUS detector
bug: a per-frame adaptive threshold that admitted clutter at low concentration, fixed by the
`'fixed'` method). This one documents what the fix actually did, why the results are
nonetheless sound, and the several red herrings that cost a full day of investigation.
**Read this before touching the detection threshold or the PSF.**

All numbers were recomputed directly from the source data by independent verification scripts
(kept in the analysis workspace, not distributed here: they hard-code machine-local acquisition
paths). Everything needed to reproduce them is in this repository plus the raw data.

---

## TL;DR

1. **The shipped detection threshold is `prctile(envPool, 99.9)`, not a calibrated knee.** The
   sweep cannot reach its own tolerance, silently takes a fallback, and returns its own search
   ceiling. All 6 (dataset x domain).
2. **The counts are still sound. `beta ~ 0.508` stands.** The Bg pedestal subtraction is
   unbiased because the noise floor is concentration-INDEPENDENT (measured). The bug costs
   PRECISION, not ACCURACY.
3. **DO NOT "fix" this by raising the threshold to a false-alarm-free point.** That introduces
   AMPLITUDE-SELECTION BIAS and fabricates a slope increase. Details in section 4. This is the
   single most important warning in this document.
4. **The PSF is >= 0.245 mm (diffraction floor) and ~0.25 mm measured.** The SUSHI kernel
   (0.2456) is correct. **0.349 and 0.252 are both retracted** (see section 6).

## 1. The sweep returns its own ceiling

`REANALYSIS_RUNNER.m / local_sweep_threshold` draws candidates as

```matlab
lo  = prctile(envPool, 50);
hi  = prctile(envPool, 99.9);      % ceiling, over the WHOLE FOV
cand = linspace(lo, hi, cfg.thrSweep.nThr);
```

`envPool` is subsampled over the whole field of view. The tube is a small fraction of that
FOV, so p99.9 of the pooled envelope sits BELOW the in-tube peak amplitudes. No candidate ever
reaches `cfg.thrSweep.tolLocPerFrame = 0.02`, so `lastDirty` is the last index, `okIdx`
overflows, and the code takes

```matlab
[~, okIdx] = min(falseRate);   % none stays clean: take the cleanest
```

which returns the LARGEST candidate, i.e. `hi` itself. **The shipped threshold is exactly
`prctile(envPool, 99.9)`.** The `WARNING: no candidate stayed <= tol` goes to a diary and was
never seen. The plot then labels the result `knee=` unconditionally.

**Evidence (6/6 combinations):**

| dataset | domain | shipped thr | measured p99.9 | Bg false rate at shipped | ratio |
|---|---|---|---|---|---|
| jun23 | PI | 420.4 | 418.8 | 18.00 /frame | 1.004 |
| jun23 | fundamental | 420.6 | 417.0 | 17.91 /frame | 1.009 |
| jun23 | singlepol | 298 | 295.3 | 18.20 /frame | 1.009 |
| apr17 | PI | 419 | 419.9 | 14.14 /frame | 0.998 |
| apr17 | fundamental | 423.4 | 421.6 | 12.40 /frame | 1.004 |
| apr17 | singlepol | 294.6 | 300.0 | 15.82 /frame | 0.982 |

Every shipped threshold is within 2% of its own p99.9. Design tolerance 0.02 false loc/frame;
delivered 12-18, ~1000x over. Corroborations: the "knee" is the last data point on every
`thr_sweep_*.png`; jun23 (420.6) and apr17 (423.40) agree to 0.7% despite thresholds being
dataset-local for sessions with different gain; and singlepol lands at ~295 rather than ~420,
i.e. a different absolute value but the same p99.9 relationship, which rules out "420 is a
physical constant".

## 2. Why the counts are sound anyway

**The noise floor is concentration-INDEPENDENT** (`noise_conc_dependence.m`, jun23
fundamental, in-tube peaks/frame in bands far below any plausible bubble):

| band | beta vs conc | rung/BGTF at L1 | rung/BGTF at U5 |
|---|---|---|---|
| [150,300) deep noise | -0.101 | 0.91x | 0.55x |
| **[300,420) noise below the shipped thr** | **-0.016** | **1.10x** | **1.04x** |
| [420,612) admitted by the shipped thr | +0.194 | 2.38x | 7.22x |
| [612,inf) above the Bg ceiling | +0.744 | 193x | 6752x |

Band 2 is flat across a **147x** concentration range. (Band 1 declines because bright real
peaks NMS-block dim noise peaks: depletion, not scaling. Bands 3 and 4 rise because that is
where real bubbles live.)

**Therefore the Bg pedestal is an unbiased estimator of rung-block false alarms**, and
subtraction removes the noise exactly. Verified at L1, jun23 fundamental:

```
raw@420.6    = 46.43/frame   (= 42.575 in [420,612) + 3.855 in [612,inf))
BGTF@420.6   = 17.93/frame   (= 17.907 + 0.020)
loc_sub      = 46.43 - 17.93                    = 28.50
real bubbles = (42.575-17.907) + (3.855-0.020)  = 28.50    <- identical
```

`loc_sub` is a correct count of real bubbles above the threshold. The design was sound; the
threshold sitting inside the noise costs precision (variance), not accuracy (bias).

## 3. Independent validations

- Shipped `Rd.pedLoc(combinedTube) = 19.06` vs an independent re-measurement of 17.93 (6%,
  explained by an 800-frame subsample vs the readout's full ~4617).
- A from-scratch re-count reproduces the published per-block loc rates at the shipped
  threshold to **median +2.9%, max 7.4%** across all 15 jun23 rungs.
- Eli reviewed 60 fps video of L5B1/M4B1 (`bmode_movie.m`) with detections overlaid and
  confirmed the localizations track real bubble signal faithfully.

## 4. DO NOT raise the threshold to a false-alarm-free point

This is the trap. Raising the threshold to where Bg gives ~0.02 false/frame (jun23
fundamental: 612.2) counts only the BRIGHT subset, and **the bright fraction is not
concentration-independent**:

- L1: 3.855 / 46.43 = **8.3%** of detections are above 612
- U5: 135.0 / 264.3 = **51.1%**

A **6x** growth, because at high concentration bubbles overlap and interfere constructively.
So a high threshold has *relaxing selectivity* as concentration rises: it inflates the top of
the ladder and steepens the slope. Measured consequence on the classical count:

```
beta @ 420.6 (pedestal-subtracted) = 0.474      (reproduces the published 0.508)
beta @ 612.2 (false-alarm-free)    = 0.879      <- LOOKS better, IS worse
bright-fraction slope = log10(51.1/8.3)/log10(147) = +0.364
0.474 + 0.364 = 0.838 ~ 0.879                   <- the entire "improvement", accounted for
```

The apparent improvement is amplitude-selection bias and nothing else. **A lower threshold
plus an unbiased pedestal subtraction beats a higher threshold with no subtraction**, because
the noise is concentration-independent and the bubble brightness distribution is not.

Corollary for the learned methods: `infer_localizers.py --target_bg_rate` defaults to 20.0,
matching the classical operating point deliberately ("-> fair comparison"). That is correct
and should stay. Re-running the nets at 0.02 would import this bias into their slopes.

## 5. The threshold's real justification

`thr = prctile(envPool, 99.9)` is a defensible operating point, but for a reason the code
never states: it is a bias/variance compromise, not a false-alarm-free cut.

- Lower thr -> catches more dim bubbles (less selection bias) but subtracts a larger pedestal
  (more variance).
- Higher thr -> less to subtract but progressively amplitude-selected (bias).

At the shipped point, L1 recovers 28.5 signal from a 46.4 raw count on a 17.9 pedestal. Noisy
but unbiased. The tolerance-based "knee" search is the wrong instrument for this metric and
should not be resurrected.

## 6. The PSF: what is true, and what is retracted

**RETRACTED: 0.349 mm.** Mine (`psf_measured_vs_model.m` -> `psf_measured_0349.mat`, both to
be deleted). It was speckle measured at the broken threshold. It reached the committee report
and must be corrected to ~0.25. It never entered the learned pipeline's training path.

**RETRACTED: 0.252 mm.** Also mine. `psf_plateau_jun23.m` shows the pooled width climbing
monotonically with amplitude and never plateauing (0.259 -> 0.270 -> 0.281 -> 0.299 mm across
bands above 612), so any single value is a band-boundary artifact.

**TRUE:**
- **Lower bound 0.2447 mm**: the diffraction limit for this aperture. Hard floor; no aperture
  beats its own diffraction limit.
- **Measured ~0.25 mm**: the dimmest confirmed-bubble band [612,700), N=3260, gives 0.248
  (baseline-subtracted) / 0.259 (pedestal-free Gaussian + free offset).
- **Known narrowing bias**: selecting on amplitude selects well-CENTRED bubbles (a bubble near
  a pixel centre is sharper and brighter than one straddling two pixels), and each stamp still
  contains ~29 dim neighbour peaks. The literature reports real systems at ~1.1x theory
  (e.g. 276 um axial / 365 um lateral measured against a stated ~110/330 um limit). Ours
  measuring 1.01-1.06x is suspiciously good.
- **Honest statement: the PSF is ~0.25 mm, consistent with the aperture-limited prediction of
  0.245 mm; true value likely 0.245-0.28.** Quote the bracket plus the physics, not a decimal.
- **SUSHI's kernel is correct**: the shipped `psf.mat` (path set by `paths.psfFile` in
  `reanalysis/reanalysis_config.m`, loaded via `local_load_psf` at REANALYSIS_RUNNER:392)
  measures **0.2456 mm lateral / 0.0694 mm axial**, i.e. the analytical model, 3.54:1
  anisotropic. **No measured-PSF rerun is warranted.** Its entire premise was the retracted
  0.349. Verify on your own install by loading `paths.psfFile` and taking the FWHM along the
  peak row and column.

**WHY A PSF MEASUREMENT CANNOT VALIDATE BUBBLE DETECTION.** A microbubble is 1-3 um against a
~245 um resolution cell: a point scatterer ~100x smaller than the cell, so its image IS the
PSF. But speckle is PSF convolved with random scatterers, so **speckle carries the SAME
diffraction-limited width**. Any width/shape/r^2 test returns the PSF by construction whether
or not real bubbles are present. **Only amplitude discriminates.** Every failed PSF attempt in
this project (0.196, 0.244, 0.259-with-1.3px-axial, 0.349, 0.252) is a rediscovery of this.

**Anisotropy matters for figures.** 0.245 lateral x 0.069 axial = **3.5:1**: single bubbles
are HORIZONTAL STREAKS in our reconstruction, not round dots. Eli's recollection of round
bubbles in B-mode is from the **native F2 display**, which is not restricted to the 64 of 256
elements the VADA research interface reads and therefore has a ~4x tighter lateral PSF
(~0.06 mm, near-isotropic against the 0.069 axial). Our reconstruction is legitimately a worse
image than the one watched during acquisition. Both observations are correct.

## 7. Red herrings purged (each cost hours)

| where | the lie | truth |
|---|---|---|
| `detect_microbubbles.m` case `'fixed'` | "calibrated on bubble-free Bg blocks so they yield ~0 detections" | yields 12-18/frame in all 6 combos; the PEDESTAL SUBTRACTION is what makes the count valid |
| `calibrate_bg_floor.m` | "The COUNT axis counts QC-tracks" | dead code the runner explicitly rejects ("the WRONG logic and is NOT carried over"); the count is the localization RATE, no tracking, no QC |
| `calibrate_bg_floor.m` docstring | `tolLocPerFrame` default 0.05 | the code says 8 |
| `local_sweep_threshold` + `thr_sweep_*.png` | "knee=" | it is the search CEILING whenever tol is unmet, which is always |
| `psf_measured_vs_model.m`, `psf_measured_0349.mat` | 0.349 mm "measured PSF" | speckle at the broken threshold; DELETE |
| rung labels | jun23 figures alias its 15 rungs as C1..C15 | collides with **apr17's real C1..C8**. Also `jun23.blockFmt='%sB%d'` (uppercase B) vs `apr17.blockFmt='%sb%d'` (lowercase b) |

## 8. What was NOT changed, and why

The shipped threshold value is preserved exactly (`prctile(envPool,99.9)`), so all published
results reproduce. The sweep's tolerance/knee machinery is documented as abandoned rather than
repaired, because repairing it (reaching 0.02) would introduce the section-4 bias. Nothing in
`readout_*.mat`, `betas_*.mat`, or the Fig 5 CSVs changes. The count metric remains the in-tube
localization RATE, pedestal-subtracted, no tracking, no QC.
