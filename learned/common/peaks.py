#!/usr/bin/env python
# SPDX-License-Identifier: MIT
"""
SHARED peak extraction: super-res density/intensity map -> localizations.
Used IDENTICALLY by sparse deconvolution, Deep-ULM, and LISTA so the only thing that varies
across the three count-axis estimators is the map producer, not the counting.

The peak-detection threshold is the analogue of the classical detector's operating point and
MUST be set by an independent, pre-registered criterion (a fixed amplitude/SNR rule or a target
false-alarm rate calibrated on a bubble-free Bg block) -- NOT tuned to reproduce the LAT-ULM
count, or the cross-algorithm invariance claim is circular (Cross_Algorithm_Analysis_Design risk #10).
"""
import numpy as np
from scipy.ndimage import maximum_filter, label, maximum_position


def extract_peaks_sr(sr_map, thr, min_sep_sr=2, exclude_border_sr=0):
    """Local maxima of an SR map above `thr`, deduped to one peak per connected cluster.
    Returns (peaks_sr [N,2] = (z_sr, x_sr) int, amps [N]).

    Uses O(N) connected-component labelling, NOT an O(N^2) greedy-NMS loop: on a near-flat map
    the flat regions tie with the max-filter and the candidate mask lights up thousands-to-millions
    of pixels, which made the old O(N^2) loop hang (the 2026-06-26 inference hang). Connected
    components collapse each flat/clustered region to its single brightest pixel and cannot blow up."""
    win = 2 * int(min_sep_sr) + 1
    mx = maximum_filter(sr_map, size=win, mode="constant")
    cand = (sr_map == mx) & (sr_map > thr)
    if exclude_border_sr > 0:
        b = exclude_border_sr
        cand[:b, :] = cand[-b:, :] = cand[:, :b] = cand[:, -b:] = False
    lbl, n = label(cand)
    if n == 0:
        return np.empty((0, 2), int), np.empty((0,), np.float32)
    pos = np.asarray(maximum_position(sr_map, lbl, index=np.arange(1, n + 1)), dtype=int).reshape(-1, 2)
    amps = sr_map[pos[:, 0], pos[:, 1]].astype(np.float32)
    return pos, amps


def peaks_to_localizations(peaks_sr, amps, sr, dx_mm, dz_mm, x0_mm, z0_mm):
    """Convert SR-grid peaks -> localizations in the MATLAB contract order [x_mm, z_mm, amp].
    SR index (z_sr, x_sr) -> native float (z_sr/sr, x_sr/sr) -> mm via the grid origin/spacing."""
    if len(peaks_sr) == 0:
        return np.empty((0, 3), np.float32)
    z_nat = peaks_sr[:, 0] / sr
    x_nat = peaks_sr[:, 1] / sr
    x_mm = x0_mm + x_nat * dx_mm
    z_mm = z0_mm + z_nat * dz_mm
    return np.column_stack([x_mm, z_mm, amps]).astype(np.float32)


def match_localizations(rec_zx, gt_zx, tol_px):
    """Greedy nearest-neighbour matching in native px. Returns (n_match, n_rec, n_gt)."""
    n_rec, n_gt = len(rec_zx), len(gt_zx)
    if n_rec == 0 or n_gt == 0:
        return 0, n_rec, n_gt
    used = np.zeros(n_gt, bool)
    n_match = 0
    for r in rec_zx:
        d2 = (gt_zx[:, 0] - r[0]) ** 2 + (gt_zx[:, 1] - r[1]) ** 2
        d2[used] = np.inf
        j = int(np.argmin(d2))
        if d2[j] <= tol_px ** 2:
            used[j] = True
            n_match += 1
    return n_match, n_rec, n_gt
