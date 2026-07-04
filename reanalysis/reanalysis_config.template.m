function paths = reanalysis_config()
%REANALYSIS_CONFIG  Machine-local paths for REANALYSIS_RUNNER.
%
%  PORTABILITY: copy this file to reanalysis_config.m (in this same folder)
%  and set the paths below for your environment. reanalysis_config.m is
%  gitignored, so machine-specific paths (data roots, output root, code/SDK
%  paths, model PSF) never enter version control; this .template.m file is the
%  committed reference. Keep the function name as reanalysis_config after
%  copying (it must match the .m filename so MATLAB can call it).
%
%  All fields below are required by REANALYSIS_RUNNER.m.

% Code / SDK locations -------------------------------------------------------
paths.repoRoot     = '<SET ME: pipeline repo root, e.g. C:\path\to\ulm-pipeline>';
paths.codePath     = '<SET ME: folder of supporting .m code to add to the MATLAB path>';
paths.vadaPath     = '<SET ME: FUJIFILM VisualSonics VEVO F2 / VADA SDK Matlab folder>';

% Model point-spread function (UHF29x) --------------------------------------
paths.psfFile      = '<SET ME: full path to psf.mat>';

% Output root (each dataset gets a subfolder created under here) -------------
paths.outputRoot   = '<SET ME: writable output root for reanalysis results>';

% Raw data roots (one per dataset) ------------------------------------------
paths.apr17RawRoot = '<SET ME: folder of apr17 .vada raw data>';
paths.jun23RawRoot = '<SET ME: folder of jun23 .vada raw data>';
end
