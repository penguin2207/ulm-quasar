# SPDX-License-Identifier: MIT
"""Machine-local paths for the learned pipeline (Python side).

PORTABILITY: copy this file to ``learned_config.py`` (at the learned/ root, next
to this template) and set the paths for your environment. ``learned_config.py``
is gitignored, so machine-specific paths never enter version control; this
template is the committed reference. Alternatively, export the environment
variables ``LEARNED_OUTPUT_ROOT`` / ``LEARNED_RAW_ROOT`` / ``LEARNED_SCRATCH_DIR``
instead of using this file. See ``learned/README.md``. The MATLAB glue reads the
analogous ``learned_config.m``.
"""

# Reanalysis output root: holds <dataset>/readout_<ds>.mat, roi_polys_<ds>.mat,
# betas_<ds>.mat (the outputRoot of reanalysis/REANALYSIS_RUNNER.m).
OUTPUT_ROOT = r"<SET ME: reanalysis output root>"

# Raw .vada data root (parent of the per-dataset raw folders). Only needed by the
# MATLAB envelope export; pure-Python inference starts from the exported envelopes.
RAW_ROOT = r"<SET ME: raw .vada data root, or leave unset if not exporting>"

# Writable scratch dir for intermediate .npy (locrate_<ds>.npy, ceiling_stats.npy).
SCRATCH_DIR = r"<SET ME: writable scratch dir>"
