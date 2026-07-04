#!/usr/bin/env python
# SPDX-License-Identifier: MIT
"""Machine-local path resolution for the learned pipeline.

No machine-specific paths are hard-coded in tracked files. Provide them EITHER by
copying ``learned/learned_config.template.py`` -> ``learned/learned_config.py``
(gitignored) and filling it in, OR by exporting the environment variables
``LEARNED_OUTPUT_ROOT`` / ``LEARNED_RAW_ROOT`` / ``LEARNED_SCRATCH_DIR``.

  OUTPUT_ROOT  reanalysis output root (holds <dataset>/readout_<ds>.mat, roi, betas)
  RAW_ROOT     raw .vada data root (used by export_svd_envelopes.m; not needed for
               pure-Python inference, which starts from the exported envelopes)
  SCRATCH_DIR  writable scratch dir for intermediate .npy (locrate_*, ceiling_stats)

See ``learned/README.md``. Importing this module never fails on missing config;
an unset root resolves to '' so config-free imports (e.g. reading the frozen
concentration axis) still work. Call :func:`require` at a path-using entry point
to fail early with a clear message.
"""
import os
import sys

_LEARNED_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))  # the learned/ dir
if _LEARNED_ROOT not in sys.path:
    sys.path.insert(0, _LEARNED_ROOT)


def _resolve(attr, env):
    val = None
    try:
        import learned_config as _cfg  # local gitignored copy at learned/learned_config.py
        val = getattr(_cfg, attr, None)
    except Exception:
        val = None
    val = val or os.environ.get(env, "")
    return "" if str(val).startswith("<SET") else str(val)


OUTPUT_ROOT = _resolve("OUTPUT_ROOT", "LEARNED_OUTPUT_ROOT")
RAW_ROOT = _resolve("RAW_ROOT", "LEARNED_RAW_ROOT")
SCRATCH_DIR = _resolve("SCRATCH_DIR", "LEARNED_SCRATCH_DIR") or os.getcwd()


def require(name):
    """Return a configured root by name, or raise a clear setup error."""
    roots = {"OUTPUT_ROOT": OUTPUT_ROOT, "RAW_ROOT": RAW_ROOT, "SCRATCH_DIR": SCRATCH_DIR}
    val = roots.get(name, "")
    if not val:
        raise RuntimeError(
            "%s is not configured. Copy learned/learned_config.template.py -> "
            "learned/learned_config.py and set it, or export LEARNED_%s. "
            "See learned/README.md." % (name, name))
    return val
