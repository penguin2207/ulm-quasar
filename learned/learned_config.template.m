% SPDX-License-Identifier: MIT
function paths = learned_config()
%LEARNED_CONFIG  Machine-local paths for the learned-localizer pipeline (learned/).
%
%  PORTABILITY: copy this file to learned_config.m (in this same folder) and set
%  the paths below for your environment. learned_config.m is gitignored, so
%  machine-specific paths never enter version control; this .template.m file is
%  the committed reference. Keep the function name learned_config after copying
%  (it must match the .m filename so MATLAB can call it).
%
%  The MATLAB glue functions (inference/export_svd_envelopes.m,
%  inference/count_external_localizers.m) take ALL paths as arguments and are
%  portable on their own; this config supplies the machine-local roots a driver
%  or runbook builds per-block paths from. The Python side reads the analogous
%  learned_config.py (see learned/README.md).

% Pipeline repo root (this ulm-pipeline checkout) ----------------------------
paths.repoRoot   = '<SET ME: pipeline repo root, e.g. C:\path\to\ulm-pipeline>';

% Reanalysis output root: holds <dataset>/pol_cache/<tag>_POL.mat,
% roi_polys_<ds>.mat and readout_<ds>.mat -- i.e. the SAME outputRoot as
% reanalysis/reanalysis_config.m (reanalysis/REANALYSIS_RUNNER.m writes it) -----
paths.outputRoot = '<SET ME: reanalysis output root (== reanalysis_config outputRoot)>';

% Writable scratch dir for learned intermediates: exported envelope .mat,
% per-block localization .mat (<method>_<tag>_locs.mat), locrate_<ds>.npy --------
paths.scratchDir = '<SET ME: writable scratch dir for learned intermediates>';

% Python interpreter used to run the learned inference (inference/infer_localizers.py)
paths.pythonExe  = '<SET ME: python executable, e.g. C:\path\to\python.exe or python>';
end
