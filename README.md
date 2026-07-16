<p align="center">
  <img src="quasar_logo.png" alt="QUASAR Logo" width="300">
</p>

<h1 align="center">LAT-ULM & QUASAR Pipelines</h1>
<p align="center">GPU-accelerated ultrasound super-resolution for VEVO F2 / VADA</p>
<p align="center"><em>Sun Lab, Northeastern University</em></p>

---

## Overview

Localization-and-Tracking Ultrasound Localization Microscopy (LAT-ULM) and
Quantitative Ultrasound Assessment via Sparse Amplitude Recovery (QUASAR) pipelines
for plane wave pulse inversion Doppler data from the VEVO F2 / VADA system.
GPU-accelerated for NVIDIA A6000 (48 GB VRAM).

**References:**
- Hingot et al., *Sci Rep* 2019;9:2456
- Heiles et al. (PALA), *Nat Biomed Eng* 2022;6(5):605-616
- Demené et al., *IEEE TMI* 2015;34(11):2271-2285 (SVD filter)
- Tang et al., *IEEE TUFFC* 2020;67(9):1738-1751 (Kalman tracking)

## Files

```
core/                              ← main pipelines + processing functions
├── LAT_ULM_pipeline.m             ← LAT-ULM main pipeline
├── QUASAR_pipeline.m              ← QUASAR main pipeline
├── beamform_planewave_gpu.m       ← GPU DAS beamformer + Hilbert on GPU
├── svd_clutter_filter_gpu.m       ← GPU SVD clutter rejection
├── svd_clutter_filter_rsvd.m      ← seeded randomized-SVD clutter filter
├── detect_microbubbles.m          ← threshold/NCC detection + NMS
├── intensity_weighted_centroid.m  ← sub-pixel localization
├── track_microbubbles.m           ← Kalman + Hungarian tracking
├── hungarian_algorithm.m          ← optimal assignment (matchpairs if available)
├── estimate_tissue_motion.m       ← phase-correlation motion estimation
├── apply_motion_correction.m      ← Fourier shift theorem correction
├── sushi_sparse_recovery.m        ← FISTA sparse recovery (QUASAR)
├── quasar_refit.m                 ← post-LASSO / relaxed-LASSO amplitude refit
└── process_single_block.m         ← per-block processing (called by pipeline)

reanalysis/                        ← UHF29x phantom concentration-response reanalysis
├── REANALYSIS_RUNNER.m            ← entry point (cfg.step 0 = ROI check / 1 = draw+preview / 2 = full)
├── fit_concentration_response.m   ← log-log slope fits + bootstrap CIs
└── reanalysis_config.template.m   ← copy to reanalysis_config.m (gitignored, machine-local paths)

preview/                           ← data exploration + QC + parameter tuning
├── VADA_explore.m                 ← RUN FIRST: inspect data, check PI quality
├── PARAMETER_TUNER.m              ← interactive single-dataset parameter tuning
├── RESULTS_VIEWER.m               ← side-by-side result comparison
├── PREVIEW_BLOCK_QC.m             ← per-block quality control
└── ...                            ← VADA organize/xml, SVD frame inspector, block diagnostics

batch/                             ← BATCH_RUN.m (multi-dataset runner), BATCH_PARAMETER_SWEEP.m (overnight sweep)
acquisition/                       ← acq_load_block_meta.m (VADA block metadata; used by reanalysis)
cuda/                              ← optional CUDA MEX beamformer (beamform_cuda.m + beamform_pw_das.cu)
figures/                           ← ANIMATE_BUBBLES.m (track animation video generator)
docs/                              ← pipeline optimization log
learned/                           ← learned count-axis localizers: Deep-ULM + LISTA (see learned/README.md)
```

## Dependencies

- **MATLAB R2019b+** with Parallel Computing Toolbox (for gpuArray)
- **Signal Processing Toolbox** (hilbert fallback on CPU)
- **Image Processing Toolbox** (normxcorr2 for NCC detection mode)
- **VADA scripts** from VisualSonics (`VsiVadaDataRead.m`, `VsiVadaConfigRead.m`, `VsiParseXml.m`)

## GPU Acceleration

| Component | GPU Method | Speedup |
|-----------|-----------|---------|
| Beamforming | Precomputed delay tables on GPU, vectorized interpolation | ~50-100x |
| Hilbert transform | FFT-based analytic signal on gpuArray | ~10-20x |
| SVD clutter filter | gpuArray-aware `svd()` on Casorati matrix | ~5-10x |
| Detection/Localization | CPU (fast enough, small data after SVD) | 1x |

**A6000 memory budget per chunk (400 frames):**
- Delay tables: 3 angles × nPixels × nChannels × 4 bytes ≈ 2-4 GB
- RF data per frame: samples × channels × 4 bytes ≈ 2-5 MB
- IQ stack: nZ × nX × 400 × 8 bytes ≈ 1-3 GB
- Total: ~5-8 GB (well within 48 GB VRAM)

Set `config.useGPU = false` to fall back to CPU processing.

## Quick Start

### 1. Explore data (single block)
```matlab
% Edit paths in VADA_explore.m, then:
>> VADA_explore
```
Check: event ordering, PI cancellation (should be < -10 dB), beamformed preview.

### 2. Run pipeline

**Single block:**
```matlab
config.mode         = 'single';
config.dataFolder   = 'C:\path\to\data';
config.baseFilename = 'Study_2026-02-04';  % filename without .vada extension
>> LAT_ULM_pipeline
```

**All blocks in a folder:**
```matlab
config.mode       = 'folder';
config.dataFolder = 'C:\path\to\data';     % scans recursively for all .vada files
% config.baseFilename is ignored in folder mode
>> LAT_ULM_pipeline
```

Folder mode discovers all `.vada` files (including subfolders), processes each
block independently (beamform → SVD → detect → localize → track), then
accumulates all tracks for final super-resolution rendering.

**Why track per-block?** The 68-second VADA inter-block gap means no bubble
can be linked across blocks. Each block yields independent tracks. The SR
density and velocity maps accumulate across all blocks.

Progress bars show block-level and chunk-level status.

## Processing Chain

```
VADA RF [samples × channels × events]
  │
  ├── Group: 6 events → 1 compound frame (3 angles × 2 PI polarities)
  ├── Pulse inversion: sum(+,-) → bubble signal
  ├── GPU beamform: precomputed delay tables, Hilbert IQ on GPU
  ├── Coherent compound across 3 angles
  ├── GPU SVD clutter filter (Casorati decomposition)
  ├── Detection: adaptive threshold + non-maximum suppression
  ├── Localization: intensity-weighted centroid (sub-pixel)
  ├── Tracking: Kalman filter + Hungarian assignment
  └── SR render: density map + velocity map (5 μm pixels)
```

## Tuning

| Problem | Fix |
|---------|-----|
| False detections | Raise `config.det.threshold` (try 8-10) |
| Missing bubbles | Lower `config.det.threshold` (try 3) |
| Tissue residual after SVD | Raise `config.svd.cutoffLow` (try 10-20) |
| Losing bubble signal in SVD | Lower `config.svd.cutoffLow` (try 2-3) |
| Short/broken tracks | Raise `config.track.maxGapFrames` (try 5) |
| Tracks jumping vessels | Lower `config.track.maxDisp_mm` |
| Noisy velocity map | Raise `config.track.minTrackLength` (try 10+) |
| GPU out of memory | Lower `config.chunkSize` (try 200) |

## Concentration-response reanalysis (`reanalysis/`)

`reanalysis/REANALYSIS_RUNNER.m` is the entry point for the two-dataset UHF29x
phantom microbubble-ladder reanalysis. It beamforms each polarity separately
(`PI = pos + neg`, `fundamental = pos - neg`, `singlepol = pos`), applies the
seeded randomized-SVD clutter filter, and reads out, per concentration, an
in-tube **localization rate** (count axis) and a **SUSHI amplitude** (FISTA
sparse recovery followed by a post-LASSO / relaxed-LASSO amplitude refit, summed
over the ROI), then fits log-log slopes versus concentration with bootstrap
confidence intervals. Set `cfg.step` to 0 (quick ROI check), 1 (draw ROIs +
preview), or 2 (full beamform + analysis).

### Configuration / portability

Machine-specific paths (data roots, output root, code/SDK paths, model PSF) are
**not** hard-coded in the runner. Copy `reanalysis/reanalysis_config.template.m`
to `reanalysis/reanalysis_config.m` and set the paths for your environment;
`reanalysis_config.m` is gitignored so local paths never enter version control.

The raw-data reader `VsiVadaDataRead` is part of the FUJIFILM VisualSonics VEVO
F2 / VADA SDK; it is **not** distributed here and must be on the MATLAB path, so
the pipeline does not run standalone without that SDK and the raw `.vada` data.
The XML helpers `read_vada_xml_params` / `get_param` are currently local
sub-functions of the pipeline files; extracting them into standalone `.m` files
would make the runner more portable, but is not required for the current
workflow.

## Learned count-axis localizers (`learned/`)

The `learned/` module provides learned localizers for the concentration **count
axis**: a Deep-ULM convolutional network and a deep-unfolded **LISTA** network,
alongside a classical **FISTA** sparse-deconvolution baseline. Each takes the
fundamental-domain reconstruction, produces a localization map, and reads out a
per-frame localization rate for cross-algorithm concentration-response
comparison. Model weights are system-specific (tied to the measured PSF) and are
**not** distributed here: retrain per system. See `learned/README.md` for
datagen, training, inference, and its own gitignored `learned_config.{m,py}`
path setup.

Method references: van Sloun et al. 2021 (Deep-ULM); Gregor & LeCun 2010 with
Luijten et al. 2023 (LISTA); Beck & Teboulle 2009 (FISTA).

## License

MIT. See [`LICENSE`](LICENSE). Copyright (c) 2026 Eli Wirth-Apley.

## Acknowledgments

Developed by Eli Wirth-Apley (Sun Lab, Northeastern University) with assistance
from Anthropic's Claude (Claude Code) for code generation, debugging, and
refactoring under the author's direction and review. Primary method references
are listed above and in the source-file headers.
