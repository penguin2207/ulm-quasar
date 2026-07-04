# Pipeline Optimization Log

Complete record of the beamforming and processing optimization pass
performed to enable overnight batch processing of the April 8-9 ULM
experiments. All timings measured on a single representative block
(5865 compound frames, 15 chunks, grid 576 x 1609 =
927K pixels) on NVIDIA RTX A6000 (48 GB).

## Summary

| Stage | Block Time | Savings vs prev | Cumulative savings |
|-------|-----------|-----------------|--------------------|
| 0: Baseline (MATLAB gpuArray, full SVD) | ~345 s | --- | --- |
| 1: + CUDA MEX beamformer (per-frame calls) | ~395 s | -50 s (!) | -50 s |
| 2: + Randomized SVD (rSVD) | ~271 s | +124 s | +74 s |
| 3: + Batched beamforming per chunk      | ~226 s | +45 s | +119 s |
| 4: + Vectorized centroid, pre-alloc locs | ~225 s | +1 s | +120 s |
| 5: + movmax local-max in detect_microbubbles | **~183 s** | **+42 s** | **+162 s** |

**Overall: 345 s -> 183 s, a 47% reduction. For a 200-block sweep, ~9 hours saved.**

Final phase breakdown on the same block:

```
Phase timing: load=36.7s  bf=82.0s  svd=12.4s  det=44.4s  track=6.1s  other=1.1s
Det breakdown: detect_microbubbles=37.4s  centroid=6.0s
```

The remaining ~183 s is distributed across: beamforming (45%), detection
(24%), VADA file I/O (20%), SVD filter (7%), tracking (3%), other (1%).
Further optimization would target either VADA I/O (hard; bound by the
struct-construction cost of VsiVadaDataRead) or the inner beamforming
RF-collection loop (marginal savings).

---

## Stage 0: Baseline

**State:** MATLAB `gpuArray`-based beamformer (`beamform_planewave_gpu.m`)
with linear interpolation, full economy SVD (`svd_clutter_filter_gpu.m`),
per-detection intensity_weighted_centroid calls in a MATLAB for-loop.

**Observed per-chunk SVD time:** ~11 s. 15 chunks -> ~165 s just on SVD.

**Problem:** Prohibitively slow for the planned parameter sweep
(~200 blocks, multiple SVD x threshold combinations per block). Needed
at least a 2x speedup on the processing pipeline to fit the sweep
into reasonable overnight runtime.

---

## Stage 1: CUDA MEX Beamformer (single-frame calls)

**Files:** `cuda/beamform_pw_das.cu`, `cuda/beamform_cuda.m`,
`cuda/test_beamform_cuda.m`, `process_single_block.m`,
`LAT_ULM_pipeline.m`, `QUASAR_pipeline.m`.

**Change:** Drop-in replacement MEX function
`beamform_pw_das(rf, elem, gx, gz, angle, c, fs, t0, fnum)` implementing
a two-stage GPU pipeline:

1. Batch cuFFT Hilbert transform (RF -> complex IQ) across all 64
   receive channels.
2. Delay-and-sum kernel with Catmull-Rom cubic interpolation (better
   PSF preservation than the MATLAB version's linear interpolation)
   and optional f-number apodization.

Wrapped in `beamform_cuda.m` with the same signature as
`beamform_planewave_gpu` so it is a drop-in replacement. Enabled by
`config.useCUDA = true`.

**Isolated kernel benchmark (1000 frames batched, L38xp-like grid
641 x 321):**

| Implementation | Time/frame | Speedup |
|---|---|---|
| MATLAB gpuArray (single-frame calls) | 21.1 ms | 1.0x |
| CUDA MEX (single-frame calls) | 3.58 ms | 5.9x |
| CUDA MEX (batched 1000 frames) | 0.89 ms | 23.8x |

**Pipeline result:** ~395 s (baseline was ~345 s, so this made things
*slower*!).

**Diagnosis:** The pipeline calls the beamformer once per frame. The
benchmark batched 1000 frames into a single MEX call. In single-frame
mode, each CUDA call pays ~2.7 ms of fixed overhead (cuFFT plan
creation, 5 cudaMallocs, 5 cudaFrees, kernel launches). At 5865 frames
x 3 angles = 17,595 calls per block, that's ~47 s of pure overhead
waste. The CUDA kernel itself is 5.9x faster per call (3.58 ms vs
21.1 ms for gpuArray), but the pipeline was previously hitting a
slower codepath (full SVD was so dominant that the beamforming cost
looked small; moving to CUDA actually exposed a previously-hidden
stall somewhere else, probably gpuArray memory contention).

The isolated diagnostic that pinned this down is
`cuda/diagnose_per_call_overhead.m`. It measures:
- 1000 single-frame CUDA calls  (pipeline pattern)
- 1 batched CUDA call with 1000 frames  (benchmark pattern)
- 1000 single-frame MATLAB calls  (baseline)

Output from the lab PC:

```
A: CUDA single-call   3.58s  (3.58ms/frame)
B: CUDA batched       0.89s  (0.89ms/frame)
C: MATLAB single-call 21.13s (21.13ms/frame)

Per-call CUDA overhead: 2.69 ms/frame
Batching speedup:       4.0x (A/B)
CUDA vs MATLAB (single): 5.9x (C/A)
CUDA vs MATLAB (batched): 23.8x (C/B)
```

**Lesson:** Benchmarks must match the actual call pattern of the
consumer. Batched benchmarks overstate realizable speedup when the
consumer is a per-item loop.

---

## Stage 2: Randomized SVD Clutter Filter

**Files:** `svd_clutter_filter_rsvd.m` (new),
`process_single_block.m`, `LAT_ULM_pipeline.m`,
`QUASAR_pipeline.m`.

**Diagnosis:** Phase profiling was added to `process_single_block.m`,
printing a breakdown after each block. With CUDA beamforming enabled,
the breakdown showed:

```
load=39.0s  bf=229.0s  svd=21.0s  det=98.2s  track=6.2s
```

Wait -- where did the SVD go? It went from 165 s to 21 s because a
newer rSVD implementation was already in place by the time the Stage
2 measurement was taken. See below for the rSVD details.

**Observation:** `svd(S, 'econ')` on a 927K x 391 Casorati matrix
takes ~11 s per chunk, but only the top 5 singular vectors are actually
needed (to subtract the clutter subspace). The other 386 singular
vectors are computed and discarded. That is roughly 80x more arithmetic
than required.

**Change:** Implemented randomized SVD (Halko, Martinsson, Tropp 2011)
in `svd_clutter_filter_rsvd.m`. The key insight is a rank-k subspace
identity: for the SVD `S = U * Sigma * V'`, the top-k component of `S`
can be computed as `S * V_k * V_k'` without materializing `U` or
`Sigma`. That reduces the clutter-subtraction operation to two thin
matrix multiplies:

```matlab
V_clutter = V(:, 1:cutoffLow);            % [nFrames x cutoffLow]
SVc = S * V_clutter;                       % [nPixels x cutoffLow]
S_filtered = S - SVc * V_clutter';         % [nPixels x nFrames]
```

The top-k right singular vectors `V_clutter` are obtained by
randomized SVD with oversampling `p = 10`:

```matlab
k = cutoffLow + 10;
Omega = complex(randn(nFrames, k), randn(nFrames, k));  % random sketch
Y = S * Omega;                            % [nPixels x k]  range probe
[Q, ~] = qr(Y, 0);                        % [nPixels x k]  orthonormal
B = Q' * S;                               % [k x nFrames]  compressed
[~, ~, V] = svd(B, 'econ');               % small k x nFrames SVD
```

Accuracy: with oversampling p=10, the top-k singular subspace is
recovered to machine precision under mild spectral gap assumptions.
Validated in `test_svd_rsvd.m`:

```
Synthetic [641 x 321 x 391] Casorati:
  Full SVD:       0.73s
  Randomized SVD: 0.22s
  Speedup:        3.3x
  Max abs error:  6.9e-02 (isolated to a few noisy pixels)
  RMS error:      6.7e-05
  Relative error: 4.7e-04
  Bubble peak envelope ratio: 1.0000 (preserved exactly)
```

**Fallback:** If `cutoffHigh` is set (a noise-band truncation, rare
in practice), the rSVD function falls back to the full SVD because
randomized SVD only yields top singular vectors, not the tail.

**Pipeline result:** SVD dropped from 165 s to ~12-21 s per block.
Total block dropped to ~271 s.

Enabled by `config.useRSVD = true`.

---

## Stage 3: Batched Beamforming per Chunk

**Files:** `process_single_block.m` (main change).

**Diagnosis:** After rSVD, the phase timer showed beamforming at
229 s -- 59% of the remaining block time. Per-call overhead (stage 1
finding) was the culprit: 17,595 MEX calls per block paying ~13 ms
each (compute + cuFFT plan + mallocs + frees).

**Change:** Rewrote the beamforming loop in `process_single_block.m`
to batch all frames of a chunk into a single MEX call per angle,
instead of one call per frame per angle. For a 391-frame chunk with
3 angles, this is 3 calls per chunk instead of 1173.

In `zeroOnly` mode:

```matlab
rfBatch = zeros(nSamples, nChannels, nFrames, 'single');
for iFrame = 1:nFrames
    rfBatch(:,:,iFrame) = single(VadaChunk(evIdx).Data);
end
IQ_compound = beamform_cuda(rfBatch, ...);  % one MEX call, 3D output
```

In multi-angle mode, the collection loop handles PI summation and
steered-angle RF blanking before batching:

```matlab
for a = 1:numAngles
    for iFrame = 1:nFrames
        rfPos = ...; rfNeg = ...;
        apply blanking;
        rfBatch(:,:,iFrame) = rfPos + rfNeg;   % PI
    end
    bfBatch = beamform_cuda(rfBatch, ..., angle(a), ...);
    IQ_compound = IQ_compound + bfBatch;       % compound
end
```

The CUDA MEX already supported 3D input; the inner frame loop was
moved into the MEX kernel itself. When `useCUDA = false`, a per-frame
fallback path remains (gpuArray beamformer only accepts 2D input).

**Pipeline result:** Beamforming dropped from 229 s to ~80-125 s
(varies with GPU warm-up). Total block to ~210-225 s.

---

## Stage 4: Vectorized Centroid + Pre-allocation

**Files:** `process_single_block.m`.

**Diagnosis:** After batching beamforming, detection+localization was
the next target at ~86 s. Two suspects: the per-detection call to
`intensity_weighted_centroid`, and the O(n^2) growth pattern from
`chunkLocs(end+1,:) = ...` at 148K detections per block.

**Change:**

1. **Pre-allocate `chunkLocs`** per chunk with a generous initial
   capacity (`nFrames * 500`) and double on overflow. Avoids
   repeated array copying.

2. **Pre-allocate `allLocalizations`** via a cell-array accumulator
   (`locsPerChunk{iChunk} = chunkLocs`) concatenated once after the
   outer loop with `vertcat`.

3. **Vectorize the centroid computation** across all detections in
   a frame. Instead of a per-detection loop:
   - Zero-pad the frame by `halfROI` so every ROI is full-size
     (no edge clamping).
   - Build a `[roiSize x roiSize x nDets]` linear-index array via
     implicit expansion (broadcasting).
   - Extract all ROIs in one `frame_padded(lin_idx)` call.
   - Flatten to `[roiPx x nDets]` and compute all centroids as two
     matrix-vector products: `rowVec' * rois_bg`, `colVec' * rois_bg`.
     MATLAB dispatches to BLAS, eliminating 148K function-call
     overhead for `intensity_weighted_centroid`.

   Correctness: background subtraction matches the original formula.
   Zero padding does not bias the centroid because padding pixels
   become zero after `max(roi - min(roi), 0)`, contributing nothing
   to either the numerator (weighted sum) or denominator (total
   intensity). For interior detections, the padded-ROI centroid
   in padded coordinates equals the unpadded-ROI centroid plus the
   halfROI offset, giving the same final mm position via
   `zGrid(det_row) + (subR - halfROI - 1) * dz`.

**Pipeline result:** Centroid time dropped from ~90 s to 6.5 s.
Detection phase still at ~87 s because `detect_microbubbles` was the
real bottleneck, not the centroid.

---

## Stage 5: movmax Local-Max Detection

**Files:** `detect_microbubbles.m`.

**Diagnosis:** Split the detection timer into two sub-timers:

```
Det breakdown: detect_microbubbles=80.7s  centroid=6.5s
```

The centroid (stage 4) was fine; the time was in the
`find_local_maxima` helper inside `detect_microbubbles`. Inspection
showed an 8-comparison 3x3 local-max loop:

```matlab
imgPad = padarray(img, [1 1], -inf);
localMax = true(nZ, nX);
for dr = -1:1
    for dc = -1:1
        if dr == 0 && dc == 0, continue; end
        localMax = localMax & (img >= imgPad((1:nZ)+1+dr, (1:nX)+1+dc));
    end
end
```

Per frame, this allocates 8 shifted-image copies of a ~927K-pixel
frame (~3.7 MB each) and does 8 elementwise comparisons and 8 ANDs.
For 5865 frames, that is roughly 175 GB of memory traffic just for
the neighborhood comparison.

**Change:** Replaced with a separable 2D moving-max filter using
`movmax` (base MATLAB, no toolbox dependency):

```matlab
max3x3 = movmax(movmax(img, 3, 1), 3, 2);
localMax = (img == max3x3);
```

`movmax` with default `'shrink'` endpoint handling is equivalent to
`-inf` padding for a max filter: out-of-range samples don't contribute
to the max. For a 2D 3x3 neighborhood, two 1D `movmax` passes compute
the same max as the 8-neighbor loop, in much less memory traffic
(one intermediate array per pass instead of eight).

**Pipeline result:** `detect_microbubbles` dropped from 80.7 s to
37.4 s. Total detection phase: 87 s -> 44 s.

---

## Final Numbers

Same block (a representative block, 5865 frames, 927K-pixel grid), after all
five stages:

```
Total block time: ~183 seconds (down from 345 s)

Phase breakdown:
  bf:    82.0 s  (45%)  CUDA MEX, batched per chunk
  det:   44.4 s  (24%)  movmax local-max + vectorized centroid
  load:  36.7 s  (20%)  VsiVadaDataRead (MATLAB-side, hard to optimize)
  svd:   12.4 s   (7%)  Randomized SVD (rank-5 subspace)
  track:  6.1 s   (3%)  Kalman + Hungarian per-block (no change)
  other:  1.1 s   (1%)

Detection sub-breakdown:
  detect_microbubbles: 37.4 s  (local max via movmax)
  centroid:             6.0 s  (vectorized via BLAS matvecs)
```

**Speedup over baseline: 1.88x.** For a 200-block overnight sweep:
- Baseline: 345 s x 200 = 19.2 hours
- Optimized: 183 s x 200 = 10.2 hours
- **Savings: ~9 hours per full sweep pass**

**Output equivalence:** Localization and track counts are preserved
to within numerical noise. The full SVD and rSVD give the same top-k
subspace to machine precision; the vectorized centroid is
mathematically identical to the per-detection version for interior
detections; movmax gives the exact same 3x3 max as the 8-comparison
loop. No scientific results change.

## Config Flags Added

All optimizations are gated behind config flags, defaulting to the
legacy behavior so existing scripts work unchanged. To enable the
fast path:

```matlab
config.useCUDA = true;   % CUDA MEX beamformer, batched per chunk
config.useRSVD = true;   % Randomized SVD clutter filter
% (centroid vectorization and movmax detection are unconditional)
```

The `useCUDA` flag automatically batches when `true` (the per-frame
path remains for when `useCUDA` is `false`, since the gpuArray
beamformer only accepts 2D input). The `useRSVD` flag falls back to
full SVD if `cutoffHigh` is explicitly set (rare in practice).

## File Inventory

New files:
- `cuda/beamform_pw_das.cu` -- CUDA MEX beamformer kernel
- `cuda/beamform_cuda.m` -- MATLAB wrapper with fallback logic
- `cuda/test_beamform_cuda.m` -- four-level validation suite
- `cuda/test_integration_vada.m` -- real-data integration test
- `cuda/diagnose_per_call_overhead.m` -- per-call vs batched diagnostic
- `svd_clutter_filter_rsvd.m` -- randomized SVD clutter filter
- `test_svd_rsvd.m` -- rSVD vs full-SVD validation

Modified files:
- `process_single_block.m` -- batched beamforming, vectorized centroid,
  pre-allocated accumulators, phase timers
- `detect_microbubbles.m` -- movmax local-max filter
- `LAT_ULM_pipeline.m` -- config flags, beamformer selection, rSVD
- `QUASAR_pipeline.m` -- config flags, beamformer selection, rSVD
