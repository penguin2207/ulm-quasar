# Learned microbubble localizers (Deep-ULM + LISTA)

Learned super-resolution localizers for the concentration-response reanalysis
count axis: a **Deep-ULM** convolutional network and a **deep-unfolded LISTA**
network, plus a classical **FISTA sparse-deconvolution** baseline. Each is a
drop-in replacement for the localization-map producer only — it consumes the same
SVD-clutter-filtered envelope the classical LAT-ULM detector sees, so the count
comparison isolates the localization algorithm ("algorithm axis").

## Layout

```
learned/
├── deepulm/model.py                 Deep-ULM CNN
├── lista/model.py                   deep-unfolded LISTA
├── deconv/fista_deconv.py           classical FISTA sparse-deconvolution baseline
├── datagen/synth.py                 synthetic training-data generator
├── psf/
│   ├── measure_psf.py               measure the system point-spread function
│   └── build_training_psf.py        assemble the training PSF
├── train.py                         train Deep-ULM / LISTA on synthetic data
├── diagnose_nets.py                 sanity / diagnostic plots
├── inference/
│   ├── datasets.py                  frozen per-dataset spec (rungs, conc axis, ROIs, paths)
│   ├── paths.py                     machine-local path resolution
│   ├── infer_localizers.py          run a net on real envelopes -> localization .mat
│   ├── count_locrate.py             in-tube localization-rate count + run_sweep.py driver
│   ├── run_sweep.py                 sweep all blocks x methods
│   ├── ceiling_analysis.py          breakpoint + per-window slopes + bootstrap 95% CIs
│   ├── export_svd_envelopes.m       MATLAB: SVD-filtered envelope frames from a polarity cache
│   └── test_reanalysis_units.py     unit checks against the classical convention
├── figures/
│   ├── fig5_house.py, fig5_pivsfund.py     cross-algorithm Fig 5 plots
│   ├── fig5_region_slopes.csv              per-region slopes + 95% CIs (block bootstrap)
│   └── fig5_combined_tube_data_decaycorr.csv   per-rung decay-corrected values + CIs
├── learned_config.template.m        machine-local paths for the MATLAB glue (copy -> learned_config.m)
└── learned_config.template.py       machine-local paths for the Python side (copy -> learned_config.py)
```

## Inputs — the fundamental reconstruction

The localizers operate on the **fundamental-domain** reconstruction produced by
`reanalysis/REANALYSIS_RUNNER.m`: the polarity-separated (`fundamental = pos - neg`),
seeded randomized-SVD clutter-filtered (cutoff 8, seed 12345) envelope `|IQf|`.
`export_svd_envelopes.m` writes those envelope frames from a block's polarity
cache, so the nets, the deconvolution baseline, and the classical detector all see
the identical filtered input. The fundamental (linear) domain is used rather than
pulse-inversion because PI's nonlinear-contrast steepening would confound a
cross-algorithm slope comparison.

## Method

- **Deep-ULM** — a CNN that maps the low-resolution envelope to a high-resolution
  localization map (van Sloun et al. 2021). Trained per system on synthetic
  microbubble fields convolved with the *measured* system PSF.
- **LISTA** — a deep-unfolded iterative sparse-coding network (Gregor & LeCun
  2010; applied to ultrasound signal processing in Luijten et al. 2023) that
  learns a fast approximation of sparse recovery for the same map.
- **FISTA deconvolution** — a classical, non-learned sparse-deconvolution baseline
  (Beck & Teboulle 2009) for comparison.
- **Count** — the in-tube localization *rate* (detections in the ROI polygon /
  nFrames), background-pedestal-subtracted, with no tracking and no QC. This is
  the identical convention to the classical reanalysis runner, so the algorithms
  are compared on the same footing.

## Model weights: not distributed (retrain per system)

Trained weights are **not** committed and **not** distributed. The networks are
trained on synthetic data built from a *specific system's* point-spread function
(UHF29x here); the weights do not transfer across probes, systems, or imaging
depths. To use the learned localizers on your own data:

1. `psf/measure_psf.py` — measure your system PSF.
2. `psf/build_training_psf.py` — assemble the training PSF.
3. `datagen/synth.py` (+ `train.py`) — generate synthetic data and train.

Interim checkpoints are not distributed with this repository (they are specific
to this system's PSF); retrain per the steps above for any other system.

## Configure (no hard-coded paths)

Machine-specific paths are not committed. Provide them either way:

- **Python:** copy `learned_config.template.py` → `learned_config.py` and fill it
  in, **or** export `LEARNED_OUTPUT_ROOT` / `LEARNED_RAW_ROOT` / `LEARNED_SCRATCH_DIR`.
- **MATLAB glue:** copy `learned_config.template.m` → `learned_config.m`.

Both `learned_config.m` and `learned_config.py` are gitignored. `OUTPUT_ROOT` is
the same reanalysis output root as `reanalysis/reanalysis_config.m`. Resolution
logic lives in `inference/paths.py`. (The two MATLAB glue functions also accept
all paths as arguments, so they are portable without the config.)

## How to run

The learned count axis sits downstream of the classical reanalysis; each script
documents its own usage in its header.

1. **Reanalysis** — run `reanalysis/REANALYSIS_RUNNER.m` to produce the per-block
   polarity caches + the classical count `readout` (needs the VADA SDK + raw
   `.vada`; see the top-level README).
2. **Export envelopes** (MATLAB) — `export_svd_envelopes(polCachePath, outMat,
   'fundamental', nFr, useGPU)` writes the SVD-filtered envelope frames.
3. **Train** (once per system) — `measure_psf.py` → `build_training_psf.py` →
   `train.py`.
4. **Infer** (Python) — `infer_localizers.py` runs a net (or the deconv baseline)
   on the exported envelopes → a per-block localization `.mat`.
5. **Count** — `count_locrate.py` / `run_sweep.py` → `locrate_<ds>.npy` (the
   in-tube localization rate, same convention as the classical anchor).
6. **Analyze** — `ceiling_analysis.py` → per-window slopes, breakpoints, and
   bootstrap CIs; `figures/` for the cross-algorithm Fig 5 plots and the
   `fig5_region_slopes.csv` per-region slope table.

## Dependencies

- **Python 3.10+**: numpy, scipy, h5py, torch (networks), matplotlib (figures).
- **MATLAB R2019b+** for the glue; it adds the repo's `core/` and `reanalysis/`
  to the path automatically. The upstream reanalysis needs the FUJIFILM
  VisualSonics VEVO F2 / VADA SDK (not distributed here).
- An NVIDIA GPU (A6000-class) is recommended for training.

## References

- van Sloun RJG, Solomon O, Bruce M, Khaing ZZ, Wijkstra H, Eldar YC, Mischi M.
  *Super-resolution ultrasound localization microscopy through deep learning.*
  IEEE Trans Med Imaging 2021;40(3):829–839. doi:10.1109/TMI.2020.3037790
- Gregor K, LeCun Y. *Learning Fast Approximations of Sparse Coding.* Proc. 27th
  Int. Conf. on Machine Learning (ICML) 2010:399–406.
- Luijten B, Chennakeshava N, Eldar YC, Mischi M, van Sloun RJG. *Ultrasound
  Signal Processing: From Models to Deep Learning.* Ultrasound Med Biol
  2023;49(3):677–698. doi:10.1016/j.ultrasmedbio.2022.11.003
- Beck A, Teboulle M. *A fast iterative shrinkage-thresholding algorithm for
  linear inverse problems.* SIAM J Imaging Sci 2009;2(1):183–202.
  doi:10.1137/080716542

## License

MIT. See the top-level [`LICENSE`](../LICENSE). Copyright (c) 2026 Eli Wirth-Apley.

## Acknowledgments

Developed by Eli Wirth-Apley (Sun Lab, Northeastern University) with assistance
from Anthropic's Claude (Claude Code) for code generation, debugging, and
refactoring under the author's direction and review. Primary method references
are listed above.
