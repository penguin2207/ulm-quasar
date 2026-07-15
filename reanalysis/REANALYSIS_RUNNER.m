%% REANALYSIS_RUNNER.m -- 2-dataset (apr17 + 6-23) cross-domain reanalysis runner
%
%  Built by COPYING the adversarially-verified single-dataset APR17_REANALYSIS.m
%  and wrapping its main body in a per-dataset loop. The apr17 RECIPE local
%  functions are preserved BYTE-IDENTICAL (polarity beamformer, per-domain
%  {PI,fundamental,single-pol}, SUSHI amplitude + in-tube localization-rate
%  readout, the threshold sweep, the ROI draw + per-block overlays, the true-Bg
%  pedestal subtraction, the log-log beta fit, the validation/report helpers).
%  ONLY the orchestration changed: the single apr17 config became
%  cfg.datasets = {apr17_profile, jun23_profile}; STEP 1 and STEP 2 are wrapped
%  in a per-dataset loop; the background blocks, manifest, concentrations, and
%  fit windows come from the ACTIVE profile.
%
%  NON-NEGOTIABLES (carried over from the verified base):
%    * NOTHING is taken from BATCH_REANALYSIS_2DATASET.m except the jun23 manifest
%      DATA (rung labels, nominal concs, block counts, bg labels). Its
%      calibrate_bg_floor / process_IQ_block / track_microbubbles /
%      filter_tracks_quality / QC-track count path is the WRONG logic and is NOT
%      carried over. The count metric remains the in-tube LOCALIZATION RATE
%      (per-frame detect_microbubbles 'fixed' -> ROI-mask count / nFrames; NO
%      tracking, NO QC) -- identical to apr17.
%    * The recipe (polarity beamform -> SVD8 seed 12345 -> per-domain
%      {PI=pos+neg, fundamental=pos-neg, single-pol=pos} -> SUSHI amp =
%      mean(|IQf|^2) -> FISTA lambda=0.10 -> quasar_refit post-LASSO -> ROI sum;
%      loc-rate; per-domain fixed threshold swept on the FLOWING-bg blocks then
%      held across rungs; true-Bg pedestal subtraction; per-block natural-depth
%      grid; log-log beta with 400-boot CI) is byte-identical for BOTH datasets.
%      ONLY the per-dataset profile differs.
%
%  ------------------------------------------------------------------------------
%  WORKFLOW (set cfg.step):
%    STEP 0 (OPTIONAL, after STEP 1; interactive OR headless) -- QUICK ROI CHECK:
%      Validate the drawn ROIs on EVERY block (all rungs + backgrounds, BOTH
%      datasets) BEFORE the long STEP 2. Beamforms only cfg.quickFrames (~60) frames
%      per block IN-MEMORY (no cache, no analysis), overlays the 5 ROI boxes on each
%      block's PI B-mode, and writes per-block PNGs + a per-dataset montage to
%      <ds.outDir>/quick_check/. Confirms ROI placement across differing FOVs.
%    STEP 1 (INTERACTIVE, on a DESKTOP with figure windows) -- DRAW ROIs, quick:
%      For EACH dataset AT THE OUTSET: beamform a single high-conc PREVIEW block
%      on-demand (a minute or two), MANUALLY draw the 5 physical-polygon ROIs, and
%      render the preview overlay (ROI boxes on the high-conc B-mode) so you validate
%      placement immediately. 10 polygons total in one sitting -> STOP. No long
%      beamform here -- you draw FIRST, the heavy work is STEP 2.
%    STEP 2 (UNATTENDED; apr17 then 6-23 END-TO-END, automatically):
%      For EACH dataset in turn: full polarity beamform (resume-safe; skips the cached
%      preview + any <label>_POL.mat) -> per-block ROI overlays for ALL rungs ->
%      per-domain analysis -> SAVE. apr17 finalizes + saves its outputs FIRST (review
%      it while 6-23 keeps going), then 6-23 continues automatically.
%      Per domain: sweep the detector threshold on the flowing-bg blocks -> one fixed
%      value; per block x domain x ROI readout (amplitude + loc-rate, NO tracking);
%      true-Bg pedestals; per-rung aggregation; log-log beta fits over the profile's
%      fit windows; density montages; soft target (apr17) / consistency (6-23);
%      dual-tube consistency.
%
%  RESUME: re-run. STEP 1 skips any block whose <label>_POL.mat already exists.
%  STEP 2 caches per-domain readouts (<ds>/readout_<domain>.mat) and skips
%  completed domains. Thresholds and pedestals are DATASET-LOCAL and NEVER shared
%  across datasets (different sessions/gain). ROIs are the one interactive
%  exception: on a STEP-1 resume you are PROMPTED, DEFAULT = RE-SELECT; a
%  re-select deletes that dataset's stale readout_*.mat.
%
%  Batch hardening (table stakes): software-OpenGL/painters renderer, driver +
%  per-dataset diaries, resume-safe namespaced timing accumulator, per-write-site
%  try/catch + disk-space guards, stack traces on error. White figure styling.
%
%  NOTE: spec section (6) names this file REANALYSIS_BOTH.m; the task instruction
%  names it REANALYSIS_RUNNER.m (followed here). BATCH_REANALYSIS_2DATASET.m is
%  DEPRECATED (wrong BGTF/QC logic); only its jun23 manifest data was reused.

close all; clear; clc;

% --- Software OpenGL + painters: avoid CEF/GPU renderer crashes in long runs. ---
set(groot, 'DefaultFigureRenderer', 'painters');
try
    opengl('save', 'software');   %#ok<OPGLO> deprecated but no replacement; harmless
catch
end

fprintf('============================================================\n');
fprintf('  REANALYSIS_RUNNER.m  (apr17 + 6-23; recipe byte-identical)\n');
fprintf('  Started: %s\n', datetime('now'));
fprintf('============================================================\n\n');

% =====================================================================
%   SECTION 2: CONFIGURATION
%   SHARED recipe constants live at the top level of cfg (one copy, both
%   datasets). DATASET COORDINATES live in the per-dataset profile structs.
% =====================================================================
cfg.step          = 1;                 % 1 = beamform+ROI+overlays (both); 2 = analysis from cache

% --- machine-local paths: loaded from reanalysis_config.m (gitignored) ---
%  Copy reanalysis_config.template.m -> reanalysis_config.m and set the paths
%  for your environment (data roots, output root, code/SDK paths, model PSF).
%  reanalysis_config.m is gitignored so machine-specific paths never enter
%  version control. See README "Configuration / portability".
if exist('reanalysis_config', 'file') ~= 2
    error('REANALYSIS_RUNNER:noConfig', ...
        ['reanalysis_config.m not found on the MATLAB path. Copy ' ...
         'reanalysis_config.template.m to reanalysis_config.m (same folder) ' ...
         'and set your local paths (it is gitignored).']);
end
paths             = reanalysis_config();
cfg.repoRoot      = paths.repoRoot;
cfg.codePath      = paths.codePath;
cfg.vadaPath      = paths.vadaPath;
cfg.psfFile       = paths.psfFile;     % SAME model PSF (same UHF29x probe); re-padded per dataset grid

% --- shared acquisition / beamform ---
cfg.modeName            = '.vada';
cfg.probePitch_fallback = 0.090;       % UHF29x safety net if XML pitch==0
cfg.maxFrames           = 0;           % 0 = all
cfg.useCUDA             = true;
cfg.useGPU              = true;
cfg.quickFrames         = 60;          % cfg.step=0 quick ROI check: frames beamformed per block (enough for SVD8 + a B-mode)

% --- shared SVD (FIXED) ---
cfg.svdCutoff = 8;
cfg.svdSeed   = 12345;

% --- shared domains ---
cfg.domains = {'PI','fundamental','singlepol'};

% --- shared SUSHI / QUASAR amplitude (FIXED) ---
cfg.lambda     = 0.10;
cfg.fista      = struct('maxIter',100,'tol',1e-4,'nonNeg',true,'useGPU',true);
cfg.refit      = struct('maxIterCG',50,'tolCG',1e-6,'useGPU',true);

% --- shared detector (FIXED method; thresh from per-domain sweep) ---
cfg.det.method      = 'fixed';
cfg.det.roiSize_px  = 5;

% --- shared threshold sweep (per domain, on each profile's flowing-bg blocks) ---
cfg.thrSweep.nThr           = 40;
cfg.thrSweep.tolLocPerFrame = 0.02;
cfg.thrSweep.span           = 'auto';
cfg.thrSweep.nFrMax         = 800;
cfg.thrFixed = struct('PI',NaN,'fundamental',NaN,'singlepol',NaN);   % NaN -> always sweep

% --- shared fit / bootstrap ---
cfg.nBoot      = 400;

% --- shared disk guards ---
cfg.minFreeGB_beamform = 80;           % jun23 caches (76 blocks, ~2x raw) are far larger than apr17's 27
cfg.minFreeGB_write    = 3;

% --- shared parent output root (each dataset gets a subfolder) ---
cfg.outputRoot = paths.outputRoot;

% =====================================================================
%   cfg.datasets = {apr17_profile, jun23_profile}   (apr17 FIRST)
% =====================================================================

% ===== apr17 profile (processed FIRST) =====
apr17.name        = 'apr17';
apr17.rawRoot     = paths.apr17RawRoot;
apr17.outDir      = fullfile(cfg.outputRoot, 'apr17');
apr17.blockFmt    = '%sb%d';                 % on-disk: C1b1.. (lowercase b; VERIFIED)
apr17.rungs       = struct( ...              % 8 rungs, 3 blocks each
    'label', {'C1','C2','C3','C4','C5','C6','C7','C8'}, ...
    'conc',  num2cell([2.9e5 4.4e5 6.6e5 9.9e5 1.5e6 2.2e6 3.3e6 5.0e6]), ...
    'nBlocks', num2cell(repmat(3,1,8)));
apr17.concPlaceholder = false;               % apr17 concs are final
apr17.bgFlow      = {'Bg1','Bg2'};           % sweep + true-Bg pedestal
apr17.bgStatic    = {'C5NoFlow'};            % static control (reported, NOT subtracted)
apr17.nTubes      = 2;
apr17.previewCands= {'C8b2','C8b1','C8b3','C7b2','C7b1','C6b2'};  % ROI-draw preview
apr17.roiFile     = fullfile(apr17.outDir, 'roi_polys_apr17.mat');
apr17.fitWindows  = {[3 7],[3 8],[1 8]};     % headline C3-C7, then C3-C8, C1-C8
apr17.headlineWin = [3 7];
apr17.target      = struct('roi','tubeL','locBeta',1.19,'ampBeta',0.94,'tol',0.25, ...
                           'window','C3-C7');               % SOFT
apr17.tubeNotes   = struct('L','LARGER lumen (legit higher signal)', ...
                           'R','smaller lumen', 'expectRatio','>1 const (size-driven)');
apr17.provenance  = {['Bg1/Bg2 flow-state UNCONFIRMED (flowing PBS vs no-flow); ' ...
                      'does not affect C3-C7 headline; confirm before quoting C1-C2.']};
apr17.notes       = 'C4 reads anomalously low in count+amp -- EXPECTED, do not "fix".';

% ===== jun23 profile (processed SECOND) =====
jun23.name        = 'jun23';
jun23.rawRoot     = paths.jun23RawRoot;
jun23.outDir      = fullfile(cfg.outputRoot, 'jun23');
jun23.blockFmt    = '%sB%d';                 % on-disk: M1B1.. (uppercase B; per BATCH manifest)
jun23.rungs       = struct( ...              % 15 rungs, variable block counts
    'label', {'L1','L2','L3','L4','L5','M1','M2','M3','M4','M5','U1','U2','U3','U4','U5'}, ...
    'conc',  num2cell([2.5e5 3.1e5 3.9e5 4.8e5 6.0e5  ...   % L1-L5 (recomputed to nominal)
                       7.0e5 1.15e6 1.9e6 3.0e6 5.0e6  ...   % M1-M5 (recomputed to nominal)
                       6.5e6 9.5e6 1.4e7 2.81e7 3.68e7]), ... % U1-U3 nominal; U4/U5 MEASURED DIRECT
    'nBlocks', num2cell([4 4 4 4 4, 6 6 6 6 6, 4 4 4 3 3]));
% ===================== CONC AXIS (RESOLVED 2026-06-26) =======================
% Note: dilutions were RECOMPUTED from each measured stock to hit the
% AGREED (nominal) concentrations -- per the dilution plan ("recompute all volumes
% from the measured value"). So L1-U3 sit at the nominal targets above; U4/U5 are the
% DIRECT Countess reads (2.81e7 / 3.68e7), which ran above their 2.05e7 / 3.0e7 nominal.
% Residual per-group Countess uncertainty (~10-20%, see as-run) shifts a group's
% absolute axis, not its internal slope; M1-M5 spans the L1-M2/M3-U3 stock boundary,
% so that adds a little scale uncertainty to the slope (central values = nominal).
jun23.concPlaceholder = false;
% Backgrounds = the FLOWING BGTF set (all 6; sweep + true-Bg pedestal):
jun23.bgFlow      = {'BGTF1','BGTF2','BGTF3','BGTF4','BGTF2B1','BGTF2B2'};
jun23.bgStatic    = {'M2SB1','M2SB2'};       % static-bubble control ONLY (NOT pedestal)
jun23.nTubes      = 2;
jun23.previewCands= {'U5B1','U5B2','U5B3','U4B1','U3B1','M5B1'};  % highest-conc first
jun23.roiFile     = fullfile(jun23.outDir, 'roi_polys_jun23.mat');
jun23.fitWindows  = {[6 10],[1 5],[11 15],[1 15]};  % headline M1-M5; ctx L,U,all
jun23.headlineWin = [6 10];                  % M1-M5 overlap = the cross-cal window
jun23.target      = [];                      % no locked target; consistency report instead
jun23.tubeNotes   = struct('L','into transducer 0.16 mL/min = CLEAN REFERENCE', ...
                           'R','away 0.08 mL/min; briefly on WITHDRAW early', ...
                           'expectRatio','~1 (same-size tubes, same conc)');
jun23.extra.measuredDirect = struct('U4',2.81e7,'U5',3.68e7);
jun23.extra.sourceStocks   = struct('L1_M2',2.42e7,'M3_U3',4.84e7,'U4_U5',7.03e7);
jun23.provenance  = {['conc = agreed nominal targets (dilutions recomputed from each ' ...
                      'measured stock per rung); U4/U5 = direct Countess reads 2.81e7/3.68e7.'], ...
                     ['RIGHT tube briefly on WITHDRAW early; per-tube cross-validate, ' ...
                      'LEFT (0.16) is the clean reference; never auto-exclude a block.'], ...
                     'BGTF backgrounds are genuinely FLOWING (rigorous vs apr17 Bg1/Bg2).'};
jun23.notes       = 'L1-L5 expected at/below floor; headline beta is M1-M5.';

cfg.datasets = {apr17, jun23};   % apr17 FIRST, jun23 SECOND
nDS = numel(cfg.datasets);

% =====================================================================
%   SECTION 3: PATH SETUP + SHADOW GUARD (shared, ONCE)
% =====================================================================
fprintf('[SETUP] Adding paths...\n');
if exist(cfg.vadaPath, 'dir'), addpath(genpath(cfg.vadaPath)); end
if exist(cfg.codePath, 'dir'), addpath(genpath(cfg.codePath)); end
% Add new core/acquisition/cuda/reanalysis AFTER codePath so the new core wins over older copies.
addpath(genpath(fullfile(cfg.repoRoot, 'core')));
addpath(genpath(fullfile(cfg.repoRoot, 'acquisition')));
addpath(genpath(fullfile(cfg.repoRoot, 'cuda')));     % beamform_cuda + MEX
addpath(fullfile(cfg.repoRoot, 'reanalysis'));

% Shadow guard: the seeded 5-arg svd_clutter_filter_rsvd must win. An old 4-arg
% copy breaks the seeded SVD.
if exist('svd_clutter_filter_rsvd','file') ~= 2 || nargin('svd_clutter_filter_rsvd') ~= 5
    error('APR17:ShadowedCode', ['svd_clutter_filter_rsvd resolves to "%s" (nargin=%d, ' ...
        'expected 5) -- an OLD copy is shadowing core/. Fix: run from a clean folder,\n' ...
        '   cd(''%s'')\nthen re-run. (run "which -all svd_clutter_filter_rsvd" to see all copies.)'], ...
        which('svd_clutter_filter_rsvd'), nargin('svd_clutter_filter_rsvd'), ...
        fullfile(cfg.repoRoot,'reanalysis'));
end

if ~exist(cfg.outputRoot, 'dir'), mkdir(cfg.outputRoot); end

% --- Driver diary + resume-safe namespaced timing accumulator (top level) ---
driverLog = fullfile(cfg.outputRoot, 'driver_log.txt');
diary(driverLog);
timingFile = fullfile(cfg.outputRoot, 'timing_report.mat');
if exist(timingFile, 'file')
    St = load(timingFile, 'timing');
    if isfield(St,'timing') && isstruct(St.timing)
        timing = St.timing;
        fprintf('  [timing] loaded existing timing_report.mat; preserving prior timings\n');
    else
        timing = struct();
    end
else
    timing = struct();
end

% --- Pre-flight disk guard (note: jun23 caches are far larger than apr17's) ---
fprintf('\n========== PRE-FLIGHT CHECKS ==========\n');
local_check_disk_space(cfg.outputRoot, cfg.minFreeGB_write, 'startup');
fprintf(['  [disk] NOTE: jun23 polarity caches (~2x raw, 76 blocks) dwarf apr17''s 27. '...
         'Confirm ample headroom on %s before STEP 1 Phase A.\n'], cfg.outputRoot);

% =====================================================================
%   SECTION 4: PREPARE EACH DATASET (manifest + dirs + series + resolve)
% =====================================================================
DS = cell(1, nDS);
for di = 1:nDS
    DS{di} = local_prepare_dataset(cfg, cfg.datasets{di});
end

% =====================================================================
%   BRANCH ON cfg.step
% =====================================================================
if cfg.step == 1
    % ---- STEP 1: draw ROIs AT THE OUTSET (preview) -> full beamform -> overlays ----
    tStep = tic;
    fprintf('\n##### STEP 1: ROI draw (outset) + polarity beamform + overlays (BOTH datasets) #####\n');

    % --- Phase 1: draw 5 ROIs for EACH dataset AT THE OUTSET ---
    %     Only a single high-conc PREVIEW block per dataset is beamformed on-demand
    %     here (quick), so ROIs are drawn BEFORE the long full beamform, not after.
    fprintf('\n========== STEP 1 PHASE 1: draw ROIs at the outset (each dataset) ==========\n');
    for di = 1:nDS
        st = DS{di};
        if ~st.available
            fprintf('  [%s] unavailable (no resolvable blocks); skip.\n', st.ds.name);
            continue;
        end
        try
            % locate, or beamform on-demand, a high-conc preview cache for the draw
            prevLabel = ''; prevCache = '';
            for pc = st.ds.previewCands
                cand = local_cache_path(st.polCacheDir, pc{1});
                if exist(cand, 'file'), prevLabel = pc{1}; prevCache = cand; break; end
            end
            if isempty(prevLabel)
                pvIdx = [];
                for pc = st.ds.previewCands
                    idx = find(strcmp({st.blocks.userLabel}, pc{1}) & ~[st.blocks.missing], 1);
                    if ~isempty(idx), pvIdx = idx; break; end
                end
                if isempty(pvIdx)
                    error('APR17:NoPreviewBlock', ['[%s] No resolvable preview block among ' ...
                        '{%s} to beamform for the ROI draw.'], st.ds.name, strjoin(st.ds.previewCands, ','));
                end
                fprintf('  [%s] preview block %s -> beamforming single block for ROI draw...\n', ...
                    st.ds.name, st.blocks(pvIdx).userLabel);
                tA = tic;
                timing = local_step1_beamform(st.cfg, st.blocks(pvIdx), st.polCacheDir, timing, timingFile);
                timing = accum_time(timing, [st.ds.name '_step1_preview'], toc(tA));
                try, save(timingFile, 'timing'); catch, end %#ok<CTCH>
                prevLabel = st.blocks(pvIdx).userLabel;
                prevCache = local_cache_path(st.polCacheDir, prevLabel);
            end
            if isempty(prevCache) || ~exist(prevCache, 'file')
                error('APR17:NoPreviewCache', '[%s] Preview cache missing after beamform attempt.', st.ds.name);
            end
            [roi, reSel] = local_select_rois(st.roiFile, prevCache, prevLabel, st.cfg, st.ds);
            DS{di}.roi = roi;
            % On RE-SELECT, invalidate this dataset's stale STEP-2 readout caches.
            if reSel
                stale = dir(fullfile(st.ds.outDir, 'readout_*.mat'));
                for s = 1:numel(stale)
                    sp = fullfile(st.ds.outDir, stale(s).name);
                    fprintf('  [ROI:%s] Re-selected: deleting stale readout cache %s\n', st.ds.name, sp);
                    try, delete(sp); catch ME, fprintf('  [ROI] WARN: could not delete %s (%s)\n', sp, ME.message); end
                end
            end
            % Render the preview-block overlay NOW so you validate ROI placement
            % immediately; the full all-rung overlays are produced per dataset in STEP 2.
            pvIdx = find(strcmp({st.blocks.userLabel}, prevLabel), 1);
            if ~isempty(pvIdx)
                local_render_overlays(st.cfg, st.blocks(pvIdx), roi, st.polCacheDir, st.overlayDir);
            end
        catch ME
            fprintf('  [ROI:%s] ERROR: %s -- skipping this dataset''s ROI + overlays.\n', st.ds.name, ME.message);
            local_print_stack(ME);
            DS{di}.roi = [];
        end
    end

    fprintf('\n##### STEP 1 COMPLETE (%.1f min) -- ROIs drawn + preview overlays #####\n', toc(tStep)/60);
    for di = 1:nDS
        st = DS{di};
        if ~st.available, continue; end
        fprintf('  [%s] ROIs: %s | preview overlay dir: %s\n', st.ds.name, st.roiFile, st.overlayDir);
    end
    fprintf(['  VALIDATE the preview overlays (ROI boxes on the high-conc B-mode) for BOTH datasets,\n' ...
             '  then set cfg.step=2 and re-run. STEP 2 then does, per dataset: full beamform +\n' ...
             '  all-rung overlays + analysis -- apr17 first (saved as it finishes), then 6-23 automatically.\n']);

elseif cfg.step == 2
    % ---- STEP 2: apr17 then jun23 SEQUENTIALLY off the caches ----
    tStep = tic;
    fprintf('\n##### STEP 2: per-domain analysis from cache (apr17 then jun23) #####\n');
    for di = 1:nDS
        st = DS{di};
        if ~st.available
            fprintf('  [%s] unavailable; skip.\n', st.ds.name);
            continue;
        end
        ds = st.ds; cfgDs = st.cfg; blocks = st.blocks;
        % Per-dataset diary for the heavy analysis; restore driver diary after.
        dsLog = fullfile(ds.outDir, sprintf('%s_reanalysis_log.txt', ds.name));
        diary off; diary(dsLog);
        fprintf('\n=============== STEP 2 [%s] ===============\n', ds.name);
        tDS = tic;
        try
            if ~exist(st.roiFile, 'file')
                error('APR17:NoROIs', ['roi_polys_%s.mat missing (%s). Run STEP 1 ' ...
                    '(cfg.step=1) interactively first.'], ds.name, st.roiFile);
            end
            Lr = load(st.roiFile, 'roi'); roi = Lr.roi;
            fprintf('  Loaded %d ROIs from %s\n', numel(roi.names), st.roiFile);

            % --- Full polarity beamform for THIS dataset (resume-safe; skips the STEP-1
            %     preview block + any already-cached block). apr17 runs first, then 6-23. ---
            tBf = tic;
            timing = local_step1_beamform(cfgDs, blocks, st.polCacheDir, timing, timingFile);
            timing = accum_time(timing, [ds.name '_step1_beamform'], toc(tBf));
            try, save(timingFile, 'timing'); catch, end %#ok<CTCH>

            % --- All-rung ROI overlay PNGs (validation record for EVERY block) ---
            local_render_overlays(cfgDs, blocks, roi, st.polCacheDir, st.overlayDir);

            gRef = local_load_any_grid(blocks, st.polCacheDir);
            psf  = local_load_psf(cfgDs.psfFile, gRef.dz, gRef.dx, cfgDs);

            [readout, betas] = local_step2_analyze(cfgDs, blocks, roi, psf, st.polCacheDir, ds);

            local_report_validation(betas, ds);          % apr17 soft target; jun23 context only
            local_dualtube_consistency(readout, betas, ds);  % per-tube agreement report
        catch ME
            fprintf('  [%s] STEP 2 ERROR: %s\n', ds.name, ME.message);
            local_print_stack(ME);
        end
        timing = accum_time(timing, [ds.name '_step2_analyze'], toc(tDS));
        try, save(timingFile, 'timing'); catch, end %#ok<CTCH>
        diary off; diary(driverLog);
        fprintf('  [%s] STEP 2 done (%.1f min) -> %s\n', ds.name, toc(tDS)/60, ds.outDir);
    end
    fprintf('\n##### STEP 2 COMPLETE (%.1f min) #####\n', toc(tStep)/60);
    fprintf('  Outputs under: %s\n', cfg.outputRoot);

elseif cfg.step == 0
    % ---- STEP 0: QUICK all-rung ROI check (frame-limited beamform; NO analysis) ----
    tStep = tic;
    fprintf('\n##### STEP 0: QUICK ROI CHECK -- %d frames/block, both datasets, no analysis #####\n', cfg.quickFrames);
    local_quick_roi_check(cfg, DS);
    fprintf('\n##### STEP 0 COMPLETE (%.1f min) #####\n', toc(tStep)/60);
    fprintf('  Per-block + montage overlays under each dataset''s quick_check/ subdir.\n');

else
    error('APR17:BadStep', 'cfg.step must be 0, 1, or 2 (got %g)', cfg.step);
end

diary off;


% #####################################################################
%   LOCAL FUNCTIONS
% #####################################################################

% ---------------------------------------------------------------------
%   ORCHESTRATION (new): per-dataset preparation
%   Builds the per-dataset cfg view (outputRoot=ds.outDir, fitWindows,
%   seriesDir), the manifest, the dirs, and resolves block basenames.
%   Returns a state struct; on a missing series the dataset is marked
%   unavailable (warn + skip) so one absent dataset never blocks the other.
% ---------------------------------------------------------------------
function st = local_prepare_dataset(cfg, ds)
    st.ds          = ds;
    st.available   = false;
    st.blocks      = struct([]);
    st.seriesDir   = '';
    st.nResolved   = 0;

    if ~exist(ds.outDir, 'dir'), mkdir(ds.outDir); end
    st.polCacheDir = fullfile(ds.outDir, 'pol_cache');
    st.overlayDir  = fullfile(ds.outDir, 'roi_overlays');
    for p = {st.polCacheDir, st.overlayDir}
        if ~exist(p{1}, 'dir'), mkdir(p{1}); end
    end
    st.roiFile = ds.roiFile;

    % Per-dataset cfg VIEW: recipe constants shared; coordinates per-profile.
    cfgDs = cfg;
    cfgDs.outputRoot = ds.outDir;
    cfgDs.fitWindows = ds.fitWindows;
    st.cfg = cfgDs;

    fprintf('\n----- dataset: %s -----\n', ds.name);
    for q = 1:numel(ds.provenance)
        fprintf('  [PROVENANCE:%s] %s\n', ds.name, ds.provenance{q});
    end
    if isfield(ds,'concPlaceholder') && ds.concPlaceholder
        fprintf('  [PROVENANCE:%s] concentrations are NOMINAL values (from config).\n', ds.name);
    end

    blocks  = local_build_blocks(ds.rungs, ds.blockFmt, ds.bgFlow, ds.bgStatic);
    nRungBlk = sum([ds.rungs.nBlocks]);
    fprintf('  [MANIFEST:%s] %d blocks (%d rung + %d bgFlow + %d bgStatic)\n', ...
        ds.name, numel(blocks), nRungBlk, numel(ds.bgFlow), numel(ds.bgStatic));

    try
        seriesDir = local_find_series(ds.rawRoot);
    catch ME
        fprintf('  [%s] WARNING: series not found under %s (%s). Dataset SKIPPED.\n', ...
            ds.name, ds.rawRoot, ME.message);
        return;
    end
    cfgDs.seriesDir = seriesDir;
    st.cfg = cfgDs;
    st.seriesDir = seriesDir;
    fprintf('  Series: %s\n', seriesDir);

    for k = 1:numel(blocks)
        try
            blocks(k).resolved = local_resolve_blockname(seriesDir, blocks(k).userLabel);
            blocks(k).missing  = false;
        catch
            blocks(k).resolved = '';
            blocks(k).missing  = true;
            fprintf('  WARNING: block %s not found\n', blocks(k).userLabel);
        end
    end
    nResolved = sum(~[blocks.missing]);
    fprintf('  Resolved %d / %d blocks\n', nResolved, numel(blocks));
    st.blocks    = blocks;
    st.nResolved = nResolved;
    st.available = nResolved > 0;
end

% ---------------------------------------------------------------------
%   Manifest -> flat block list  (signature generalized: rungs struct array
%   + on-disk blockFmt; apr17 '%sb%d' / jun23 '%sB%d'; tagging unchanged)
% ---------------------------------------------------------------------
function blocks = local_build_blocks(rungs, blockFmt, bgFlow, bgStatic)
    blocks = struct('userLabel',{},'rung',{},'blk',{},'kind',{}, ...
                    'conc',{},'rungIdx',{},'resolved',{},'missing',{});
    for iR = 1:numel(rungs)
        for b = 1:rungs(iR).nBlocks
            blocks(end+1) = struct( ...
                'userLabel', sprintf(blockFmt, rungs(iR).label, b), 'rung', rungs(iR).label, ...
                'blk', b, 'kind', 'rung', 'conc', rungs(iR).conc, 'rungIdx', iR, ...
                'resolved', '', 'missing', true); %#ok<AGROW>
        end
    end
    for b = 1:numel(bgFlow)
        blocks(end+1) = struct('userLabel', bgFlow{b}, 'rung', bgFlow{b}, ...
            'blk', 1, 'kind', 'bgFlow', 'conc', NaN, 'rungIdx', NaN, ...
            'resolved', '', 'missing', true); %#ok<AGROW>
    end
    for b = 1:numel(bgStatic)
        blocks(end+1) = struct('userLabel', bgStatic{b}, 'rung', bgStatic{b}, ...
            'blk', 1, 'kind', 'bgStatic', 'conc', NaN, 'rungIdx', NaN, ...
            'resolved', '', 'missing', true); %#ok<AGROW>
    end
end

function p = local_cache_path(polCacheDir, userLabel)
    p = fullfile(polCacheDir, [userLabel '_POL.mat']);
end

% ---------------------------------------------------------------------
%   Per-block NATURAL grid (common lateral aperture; per-block axial)
% ---------------------------------------------------------------------
function g = local_block_grid(meta)
    dz = meta.lambda_mm/2;  dx = meta.lambda_mm/2;
    rxElems = meta.elemPos_mm(meta.anglePairs(meta.zeroAngleIdx).rxElements);
    xMargin = 0.5;
    g.xGrid = (min(rxElems)-xMargin):dx:(max(rxElems)+xMargin);
    zTop = meta.depthOffset_mm;
    zBot = meta.depthOffset_mm + (meta.nSamples/(meta.fs_MHz*2))*(meta.c*1e-3);
    g.zGrid = zTop:dz:zBot;
    g.dx = dx; g.dz = dz;
    g.nX = numel(g.xGrid); g.nZ = numel(g.zGrid);
    g.minSep_mm = round(meta.lambda_mm, 3);     % minSep = lambda
end

% ---------------------------------------------------------------------
%   Beamformer selection
% ---------------------------------------------------------------------
function bf = local_pick_beamformer(cfg)
    if cfg.useCUDA && (exist('beamform_cuda','file')==2 || exist('beamform_cuda','file')==3)
        bf = @beamform_cuda;
    else
        bf = @beamform_planewave_gpu;
    end
end

% Robust 3-D-batch beamform: prefer one 3-D call (CUDA MEX path); if the
% beamformer only yields 2-D (gpuArray fallback), loop per frame.
function img = local_bf_apply(bf_fn, rf, rxPos, ang, xGrid, zGrid, fs, c, depthOff, dt)
    nF = size(rf, 3);
    if nF == 1
        img = single(bf_fn(rf, rxPos, ang, xGrid, zGrid, fs, c, depthOff, dt));
        return;
    end
    if strcmp(func2str(bf_fn), 'beamform_cuda')
        out = single(bf_fn(rf, rxPos, ang, xGrid, zGrid, fs, c, depthOff, dt));
        if size(out, 3) == nF, img = out; return; end
    end
    nZ = numel(zGrid); nX = numel(xGrid);
    img = complex(zeros(nZ, nX, nF, 'single'));
    for iF = 1:nF
        img(:,:,iF) = single(bf_fn(rf(:,:,iF), rxPos, ang, xGrid, zGrid, fs, c, depthOff, dt));
    end
end

% ---------------------------------------------------------------------
%   Polarity-preserving 5-angle compound beamform (THE core change)
%   Beamform pos and neg SEPARATELY per angle -> two complex stacks.
%   IQ_pos + IQ_neg == old PI-only cache (DAS + Hilbert are LINEAR), which is
%   what keeps the PI validation target valid. Domains are derived in STEP 2.
% ---------------------------------------------------------------------
function [IQ_pos, IQ_neg, g, bm] = local_beamform_block_pol(cfg, basename, frameLimit)
    if nargin < 3, frameLimit = []; end   % [] = ALL frames (STEP 1/2 unchanged); set only for the cfg.step=0 quick ROI check
    blockMeta = acq_load_block_meta(cfg.seriesDir, basename, cfg.modeName, [], cfg.probePitch_fallback);
    g = local_block_grid(blockMeta);

    nFrames = blockMeta.numCompoundFrames;
    if cfg.maxFrames > 0, nFrames = min(nFrames, cfg.maxFrames); end
    if ~isempty(frameLimit), nFrames = min(nFrames, frameLimit); end   % quick ROI check: only a few frames (caps event reads + loops)
    nEvt = nFrames * blockMeta.eventsPerFrame;

    bf_fn = local_pick_beamformer(cfg);
    VadaData = VsiVadaDataRead(cfg.seriesDir, basename, 1:nEvt, cfg.modeName);
    nS  = size(VadaData(1).Data, 1);
    nRx = size(VadaData(1).Data, 2);

    % NOTE: peak RAM is ~2x the old PI-only path -- two full complex-single stacks
    % (IQ_pos + IQ_neg) are held at once. Fine on the lab PC for apr17; watch it for
    % the heavier 6-23 sibling (more frames/rungs).
    IQ_pos = complex(zeros(g.nZ, g.nX, nFrames, 'single'));
    haveNeg = blockMeta.hasPI;
    if haveNeg
        IQ_neg = complex(zeros(g.nZ, g.nX, nFrames, 'single'));
    else
        IQ_neg = [];
        fprintf('    WARNING: block %s reports hasPI=false; single-pol only (no neg)\n', basename);
    end

    nAngPos = 0; nAngNeg = 0;            % angle-set guard (PI validity)
    for a = 1:blockMeta.numAngles
        ap    = blockMeta.anglePairs(a);
        rxPos = blockMeta.elemPos_mm(ap.rxElements);
        dt = bf_fn([], rxPos, ap.angle, g.xGrid, g.zGrid, ...
                   blockMeta.fs_MHz, blockMeta.c, blockMeta.depthOffset_mm, []);

        % --- POS polarity ---
        rfP = zeros(nS, nRx, nFrames, 'single');
        for iF = 1:nFrames
            rfP(:,:,iF) = single(VadaData((iF-1)*blockMeta.eventsPerFrame + ap.posIdx).Data);
        end
        IQ_pos = IQ_pos + local_bf_apply(bf_fn, rfP, rxPos, ap.angle, g.xGrid, g.zGrid, ...
            blockMeta.fs_MHz, blockMeta.c, blockMeta.depthOffset_mm, dt);
        nAngPos = nAngPos + 1;
        clear rfP;

        % --- NEG polarity (only if a PI pair exists) ---
        if haveNeg && ~isempty(ap.negIdx)
            rfN = zeros(nS, nRx, nFrames, 'single');
            for iF = 1:nFrames
                rfN(:,:,iF) = single(VadaData((iF-1)*blockMeta.eventsPerFrame + ap.negIdx).Data);
            end
            IQ_neg = IQ_neg + local_bf_apply(bf_fn, rfN, rxPos, ap.angle, g.xGrid, g.zGrid, ...
                blockMeta.fs_MHz, blockMeta.c, blockMeta.depthOffset_mm, dt);
            nAngNeg = nAngNeg + 1;
            clear rfN;
        end
    end
    clear VadaData;

    % Angle-set guard: PI = IQ_pos + IQ_neg is only valid (full fundamental
    % cancellation) if pos and neg were compounded over the SAME angles. apr17
    % 5-angle PI blocks all carry both; this is a latent safety guard.
    if haveNeg
        assert(nAngNeg == nAngPos, 'APR17:AngleMismatch', ...
            ['PI invalid for block %s: pos over %d angles, neg over %d -- ' ...
             'pos/neg must span the same angle set.'], basename, nAngPos, nAngNeg);
    end

    % Compact blockMeta subset for the cache.
    bm = struct('probeName',blockMeta.probeName, 'nElements',blockMeta.nElements, ...
        'pitch_mm',blockMeta.pitch_mm, 'fs_MHz',blockMeta.fs_MHz, 'c',blockMeta.c, ...
        'txFreq_MHz',blockMeta.txFreq_MHz, 'lambda_mm',blockMeta.lambda_mm, ...
        'depthOffset_mm',blockMeta.depthOffset_mm, 'nSamples',blockMeta.nSamples, ...
        'frameRate_Hz',blockMeta.frameRate_Hz, 'eventsPerFrame',blockMeta.eventsPerFrame, ...
        'numAngles',blockMeta.numAngles, 'hasPI',blockMeta.hasPI);
end

% ---------------------------------------------------------------------
%   STEP 1: beamform polarity caches (resume-safe per block)  [VERBATIM]
% ---------------------------------------------------------------------
function timing = local_step1_beamform(cfg, blocks, polCacheDir, timing, timingFile)
    fprintf('\n  ===== beamform -> pol_cache (resume-safe) =====\n');
    tA = tic;
    proc = find(~[blocks.missing]);
    for kk = 1:numel(proc)
        k = proc(kk);
        cachePath = local_cache_path(polCacheDir, blocks(k).userLabel);
        if exist(cachePath, 'file')
            fprintf('  [%d/%d] %s -- cached, skip\n', kk, numel(proc), blocks(k).userLabel);
            continue;
        end
        fprintf('  [%d/%d] %s (%s) -- polarity beamforming...\n', kk, numel(proc), ...
            blocks(k).userLabel, blocks(k).resolved);
        tBlk = tic;
        try
            [IQ_pos, IQ_neg, g, bm] = local_beamform_block_pol(cfg, blocks(k).resolved);
            local_check_disk_space(polCacheDir, cfg.minFreeGB_beamform, ...
                sprintf('POL save %s', blocks(k).userLabel));
            S = struct('IQ_pos',IQ_pos, 'IQ_neg',IQ_neg, 'nFrames',size(IQ_pos,3), ...
                       'g',g, 'blockMeta',bm);
            save(cachePath, '-struct', 'S', '-v7.3');
            fprintf('    %d frames, %.1f sec, %.1f GB\n', size(IQ_pos,3), toc(tBlk), ...
                dir(cachePath).bytes/1e9);
            clear IQ_pos IQ_neg S;
        catch ME
            fprintf('    ERROR: %s\n', ME.message);
            local_print_stack(ME);
        end
    end
    timing = accum_time(timing, 'step1_beamform', toc(tA));
    try, save(timingFile, 'timing'); catch, end %#ok<CTCH>
end

% ---------------------------------------------------------------------
%   PI log-power B-mode of a block from its polarity cache (SVD8)  [VERBATIM]
% ---------------------------------------------------------------------
function [bmode, g] = local_preview_bmode_pol(cachePath, cfg)
    if ~exist(cachePath, 'file')
        error('APR17:NoCache', 'Polarity cache missing: %s', cachePath);
    end
    L = load(cachePath, 'IQ_pos', 'IQ_neg', 'g');
    if isempty(L.IQ_neg), D = L.IQ_pos; else, D = L.IQ_pos + L.IQ_neg; end
    IQf = svd_clutter_filter_rsvd(D, cfg.svdCutoff, [], cfg.useGPU, cfg.svdSeed);
    bmode = 10*log10(mean(abs(IQf).^2, 3) + eps);
    g = L.g;
end

% ---------------------------------------------------------------------
%   MANUAL ROI draw (5 physical polygons; never hardcoded / never auto)
%   [VERBATIM draw logic; previewCands + roiFile come from ds; prompt names ds]
% ---------------------------------------------------------------------
function [roi, reSelected] = local_select_rois(roiFile, prevCache, prevLabel, cfg, ds)
    reSelected  = false;
    haveROIs    = exist(roiFile, 'file');
    interactive = usejava('desktop') && feature('ShowFigureWindows');

    if ~interactive
        % STEP 1 requires drawing; with no desktop we cannot draw.
        error('APR17:Headless', ['STEP 1 ROI draw needs a DESKTOP with figure windows. ' ...
            'Run the ROI step interactively first, then run STEP 2 headless.']);
    end

    fprintf('\n  [ROI:%s] selecting 5 physical-polygon ROIs (preview %s)\n', ds.name, prevLabel);

    % Resume: PROMPT reuse-vs-reselect, DEFAULT = RE-SELECT (never silently reuse).
    if haveROIs
        fprintf('  [ROI:%s] Saved ROI set exists: %s\n', ds.name, roiFile);
        ans_  = input('        Reuse saved ROIs? (NEVER silently reused) [y/N]: ', 's');
        if strcmpi(strtrim(ans_), 'y')
            L = load(roiFile, 'roi'); roi = L.roi;
            fprintf('  [ROI] Reusing saved ROIs.\n');
            return;
        end
        fprintf('  [ROI] Re-selecting (default).\n');
    end
    reSelected = true;

    [bmode, g] = local_preview_bmode_pol(prevCache, cfg);
    fprintf('  [ROI] Preview block: %s\n', prevLabel);

    names   = {'full','combinedTube','tubeL','tubeR','background'};
    prompts = {'FULL FOV (whole image; no draw)', ...
               'COMBINED tube (BOTH lumens) -- the PRIMARY signal ROI', ...
               'LEFT tube only (the LARGER lumen; consistency/resolvability check)', ...
               'RIGHT tube only', ...
               'BACKGROUND (drawn bubble-free region; SEPARATE from Bg1/Bg2 blocks)'};

    roi = struct('names', {names}, 'poly', {cell(1,numel(names))}, ...
                 'previewBlock', prevLabel, 'previewGrid', g);

    for i = 1:numel(names)
        fig = figure('Name', sprintf('ROI %d/%d: %s', i, numel(names), names{i}), ...
                     'Color','w', 'Position', [80 80 900 650]);
        ax = axes('Parent', fig); set(ax, 'Color','w');
        imagesc(ax, g.xGrid, g.zGrid, bmode);
        axis(ax,'image'); set(ax,'YDir','normal'); colormap(ax,hot); colorbar(ax);
        xlabel(ax,'Lateral (mm)'); ylabel(ax,'Axial (mm)');
        title(ax, sprintf('Draw ROI %d/%d: %s  (double-click to close polygon)', ...
            i, numel(names), prompts{i}), 'Interpreter','none');
        fprintf('  [ROI %d/%d] Draw: %s\n', i, numel(names), prompts{i});

        if strcmp(names{i}, 'full')
            % SENTINEL: 'full' = whole FOV, recomputed PER BLOCK from each block's
            % own grid (blocks have different natural depths). poly stays empty.
            poly = [];
            uiwait(msgbox(['FULL-FOV ROI = entire image of EACH block (no draw; ' ...
                'computed per block at mask time). Click OK.'], 'Full FOV', 'modal'));
        else
            h = drawpolygon('Color','c');
            poly = h.Position;                   % [x_mm z_mm] PHYSICAL vertices
        end
        roi.poly{i} = poly;
        close(fig);
    end

    try
        local_check_disk_space(fileparts(roiFile), cfg.minFreeGB_write, 'roi save');
        save(roiFile, 'roi');
        fprintf('  [ROI] Saved %d physical polygons -> %s\n', numel(names), roiFile);
    catch ME
        fprintf('  [ROI] WARN: could not save ROIs (%s)\n', ME.message);
    end
end

% ---------------------------------------------------------------------
%   physical-coordinate polygon -> logical mask on (xGrid,zGrid)  [VERBATIM]
% ---------------------------------------------------------------------
function mask = local_poly2mask(poly, xGrid, zGrid)
    [XX, ZZ] = meshgrid(xGrid, zGrid);
    mask = inpolygon(XX, ZZ, poly(:,1), poly(:,2));
end

% Whole-FOV rectangle (physical mm vertices) for a given block grid.
function p = local_full_poly(g)
    p = [g.xGrid([1 end end 1])', g.zGrid([1 1 end end])'];
end

% ROI rr -> logical mask on grid g. 'full' is a PER-BLOCK sentinel (true over THIS
% block's whole grid); the 4 drawn ROIs re-rasterize their physical polygon.
function mask = local_roi_mask(roi, rr, g)
    if strcmp(roi.names{rr}, 'full')
        mask = true(g.nZ, g.nX);
    else
        mask = local_poly2mask(roi.poly{rr}, g.xGrid, g.zGrid);
    end
end

% ---------------------------------------------------------------------
%   Per-block ROI overlay PNGs for ALL blocks (the validation gate)  [VERBATIM]
% ---------------------------------------------------------------------
function local_render_overlays(cfg, blocks, roi, polCacheDir, overlayDir)
    fprintf('\n  ===== rendering per-block ROI overlays =====\n');
    cols = struct('full',[1 1 1], 'combinedTube',[0 1 1], 'tubeL',[1 0 1], ...
                  'tubeR',[1 1 0], 'background',[0 1 0]);
    for k = 1:numel(blocks)
        if blocks(k).missing, continue; end
        cachePath = local_cache_path(polCacheDir, blocks(k).userLabel);
        if ~exist(cachePath, 'file')
            fprintf('  [overlay] %s: no cache, skip\n', blocks(k).userLabel);
            continue;
        end
        try
            [bmode, g] = local_preview_bmode_pol(cachePath, cfg);
            fig = figure('Visible','off', 'Color','w', 'Position',[60 60 980 720]);
            ax  = axes('Parent', fig); set(ax,'Color','w'); hold(ax,'on');
            imagesc(ax, g.xGrid, g.zGrid, bmode);
            axis(ax,'image'); set(ax,'YDir','normal'); colormap(ax,hot);
            cb = colorbar(ax); cb.Label.String = 'dB';
            for r = 1:numel(roi.names)
                if strcmp(roi.names{r}, 'full')
                    P = local_full_poly(g);          % per-block FOV box
                else
                    P = roi.poly{r};
                end
                if isempty(P), continue; end
                cc = cols.(matlab.lang.makeValidName(roi.names{r}));
                plot(ax, [P(:,1); P(1,1)], [P(:,2); P(1,2)], '-', ...
                    'Color', cc, 'LineWidth', 1.4);
            end
            if isnan(blocks(k).conc)
                ttl = sprintf('%s  (%s)  [grid %dx%d, z=%.1f-%.1f mm]', blocks(k).userLabel, ...
                    blocks(k).kind, g.nZ, g.nX, g.zGrid(1), g.zGrid(end));
            else
                ttl = sprintf('%s  %.2g MB/mL  [grid %dx%d, z=%.1f-%.1f mm]', blocks(k).userLabel, ...
                    blocks(k).conc, g.nZ, g.nX, g.zGrid(1), g.zGrid(end));
            end
            title(ax, ttl, 'Interpreter','none');
            xlabel(ax,'Lateral (mm)'); ylabel(ax,'Axial (mm)');
            outP = fullfile(overlayDir, [blocks(k).userLabel '.png']);
            local_check_disk_space(overlayDir, cfg.minFreeGB_write, ['overlay ' blocks(k).userLabel]);
            exportgraphics(fig, outP, 'Resolution', 150, 'BackgroundColor','white');
            close(fig);
            fprintf('  [overlay] %s -> %s\n', blocks(k).userLabel, outP);
        catch ME
            fprintf('  [overlay] %s ERROR: %s\n', blocks(k).userLabel, ME.message);
            local_print_stack(ME);
        end
    end
    fprintf('  Overlays written. Cyan=combinedTube, magenta=tubeL, yellow=tubeR, green=background, white=full.\n');
end

% ---------------------------------------------------------------------
%   QUICK all-rung ROI check (cfg.step==0): frame-limited beamform of EVERY
%   block (in-memory, NO cache, NO analysis) -> PI B-mode + 5 ROI-box overlay
%   PNG per block + a per-dataset montage. Lets you confirm ROI placement on
%   all rungs (incl. differing FOVs) BEFORE the long STEP 2. Reuses the verified
%   beamformer via local_beamform_block_pol(..., cfg.quickFrames).
% ---------------------------------------------------------------------
function local_quick_roi_check(cfg, DS)
    cols = struct('full',[1 1 1], 'combinedTube',[0 1 1], 'tubeL',[1 0 1], ...
                  'tubeR',[1 1 0], 'background',[0 1 0]);
    for di = 1:numel(DS)
        st = DS{di};
        if ~st.available
            fprintf('  [%s] unavailable (no resolvable blocks); skip quick check.\n', st.ds.name);
            continue;
        end
        ds = st.ds;
        if ~exist(st.roiFile, 'file')
            fprintf('  [%s] no ROI file (%s) -- run STEP 1 (cfg.step=1) first; skipping.\n', ds.name, st.roiFile);
            continue;
        end
        Lr = load(st.roiFile, 'roi'); roi = Lr.roi;
        qcDir = fullfile(ds.outDir, 'quick_check');
        if ~exist(qcDir, 'dir'), mkdir(qcDir); end
        proc = find(~[st.blocks.missing]);
        fprintf('\n=========== QUICK ROI CHECK [%s] : %d ROIs, %d blocks, %d frames/block ===========\n', ...
            ds.name, numel(roi.names), numel(proc), cfg.quickFrames);

        tiles = cell(1, numel(proc)); tlabs = cell(1, numel(proc)); tgs = cell(1, numel(proc));
        for kk = 1:numel(proc)
            k = proc(kk);
            try
                tBlk = tic;
                [IQ_pos, IQ_neg, g] = local_beamform_block_pol(st.cfg, st.blocks(k).resolved, cfg.quickFrames);
                if isempty(IQ_neg), D = IQ_pos; else, D = IQ_pos + IQ_neg; end
                IQf   = svd_clutter_filter_rsvd(D, cfg.svdCutoff, [], cfg.useGPU, cfg.svdSeed);
                bmode = 10*log10(mean(abs(IQf).^2, 3) + eps);
                clear IQ_pos IQ_neg D IQf;
                tiles{kk} = bmode; tlabs{kk} = st.blocks(k).userLabel; tgs{kk} = g;

                fig = figure('Visible','off', 'Color','w', 'Position',[60 60 980 720]);
                ax  = axes('Parent', fig); set(ax,'Color','w'); hold(ax,'on');
                local_draw_overlay(ax, bmode, g, roi, cols);
                cb = colorbar(ax); cb.Label.String = 'dB';
                if isnan(st.blocks(k).conc)
                    ttl = sprintf('%s  (%s)  [grid %dx%d, z=%.1f-%.1f mm]', st.blocks(k).userLabel, ...
                        st.blocks(k).kind, g.nZ, g.nX, g.zGrid(1), g.zGrid(end));
                else
                    ttl = sprintf('%s  %.2g MB/mL  [grid %dx%d, z=%.1f-%.1f mm]', st.blocks(k).userLabel, ...
                        st.blocks(k).conc, g.nZ, g.nX, g.zGrid(1), g.zGrid(end));
                end
                title(ax, ttl, 'Interpreter','none');
                xlabel(ax,'Lateral (mm)'); ylabel(ax,'Axial (mm)');
                outP = fullfile(qcDir, [st.blocks(k).userLabel '_quickroi.png']);
                local_check_disk_space(qcDir, cfg.minFreeGB_write, ['quickroi ' st.blocks(k).userLabel]);
                exportgraphics(fig, outP, 'Resolution', 150, 'BackgroundColor','white');
                close(fig);
                fprintf('  [%d/%d] %s -> %s (%.1fs)\n', kk, numel(proc), st.blocks(k).userLabel, outP, toc(tBlk));
            catch ME
                fprintf('  [%d/%d] %s QUICK ERROR: %s\n', kk, numel(proc), st.blocks(k).userLabel, ME.message);
                local_print_stack(ME);
            end
        end

        % --- per-dataset montage (every block in one scannable image) ---
        if any(~cellfun(@isempty, tiles))
            nT = numel(proc); gridN = ceil(sqrt(nT));
            fig = figure('Visible','off', 'Color','w', 'Position',[20 20 1700 1100]);
            tl  = tiledlayout(fig, gridN, gridN, 'Padding','compact', 'TileSpacing','compact');
            for kk = 1:nT
                ax = nexttile(tl); set(ax,'Color','w');
                if isempty(tiles{kk})
                    axis(ax,'off'); title(ax,[tlabs{kk} ' (failed)'], 'FontSize',7, 'Interpreter','none'); continue;
                end
                hold(ax,'on');
                local_draw_overlay(ax, tiles{kk}, tgs{kk}, roi, cols);
                set(ax,'FontSize',7);
                title(ax, tlabs{kk}, 'FontSize',8, 'Interpreter','none');
            end
            title(tl, sprintf(['%s quick ROI check (%d fr/block; cyan=combinedTube, magenta=tubeL, ' ...
                'yellow=tubeR, green=bg, white=full)'], ds.name, cfg.quickFrames), 'FontSize',11, 'Interpreter','none');
            outM = fullfile(qcDir, sprintf('%s_quickroi_montage.png', ds.name));
            local_check_disk_space(qcDir, cfg.minFreeGB_write, ['quickroi montage ' ds.name]);
            exportgraphics(fig, outM, 'Resolution', 150, 'BackgroundColor','white');
            close(fig);
            fprintf('  [%s] montage -> %s\n', ds.name, outM);
        end
        fprintf('  [%s] quick check done -> %s\n', ds.name, qcDir);
    end
end

% ROI-overlay drawing onto a prepared (hold-on) axes: B-mode + the 5 ROI boxes,
% matching local_render_overlays' colors + rasterization. Used by the quick check.
function local_draw_overlay(ax, bmode, g, roi, cols)
    imagesc(ax, g.xGrid, g.zGrid, bmode);
    axis(ax,'image'); set(ax,'YDir','normal'); colormap(ax,hot);
    for r = 1:numel(roi.names)
        if strcmp(roi.names{r}, 'full')
            P = local_full_poly(g);          % per-block FOV box
        else
            P = roi.poly{r};
        end
        if isempty(P), continue; end
        cc = cols.(matlab.lang.makeValidName(roi.names{r}));
        plot(ax, [P(:,1); P(1,1)], [P(:,2); P(1,2)], '-', 'Color', cc, 'LineWidth', 1.4);
    end
end

% ---------------------------------------------------------------------
%   STEP 2 helpers  [VERBATIM]
% ---------------------------------------------------------------------
function g = local_load_any_grid(blocks, polCacheDir)
    for k = 1:numel(blocks)
        if blocks(k).missing, continue; end
        cp = local_cache_path(polCacheDir, blocks(k).userLabel);
        if exist(cp, 'file')
            L = load(cp, 'g'); g = L.g; return;
        end
    end
    error('APR17:NoCacheForGrid', 'No polarity cache found to read a reference grid.');
end

% PSF: load + validate + extract a compact odd kernel; write psf_validation.png.
function psf = local_load_psf(psfFile, dz_mm, dx_mm, cfg)
    if ~exist(psfFile, 'file')
        error('APR17:NoPSF', 'PSF file not found: %s', psfFile);
    end
    P = load(psfFile); fn = fieldnames(P);
    psf0 = []; fld = '';
    for i = 1:numel(fn)
        v = P.(fn{i});
        if isnumeric(v) && ismatrix(v) && isreal(v) && all(size(v) >= 3)
            psf0 = double(v); fld = fn{i}; break;
        end
    end
    if isempty(psf0)
        % Fallback: first 2-D numeric field (even if complex -> take real part).
        for i = 1:numel(fn)
            v = P.(fn{i});
            if isnumeric(v) && ismatrix(v) && all(size(v) >= 3)
                psf0 = double(real(v)); fld = fn{i}; break;
            end
        end
    end
    if isempty(psf0)
        error('APR17:BadPSF', 'No 2-D numeric matrix field in %s', psfFile);
    end

    [Ny, Nx] = size(psf0);
    [pk, idx] = max(psf0(:));
    [pr, pc]  = ind2sub([Ny Nx], idx);

    % --- FWHM along the peak row/col ---
    fwhmLat_px = local_fwhm_px(psf0(pr, :), pc);
    fwhmAx_px  = local_fwhm_px(psf0(:, pc), pr);
    fwhmLat_mm = fwhmLat_px * dx_mm;
    fwhmAx_mm  = fwhmAx_px  * dz_mm;
    integ      = sum(psf0(:));

    fprintf('\n  ===== PSF validation =====\n');
    fprintf('  field "%s" size %s ; peak=%.4g at (row %d/%d, col %d/%d)\n', ...
        fld, mat2str([Ny Nx]), pk, pr, Ny, pc, Nx);
    fprintf('  integral=%.4g ; FWHM axial=%.3f px (%.3f mm), lateral=%.3f px (%.3f mm)\n', ...
        integ, fwhmAx_px, fwhmAx_mm, fwhmLat_px, fwhmLat_mm);

    % --- Validation figure (saved BEFORE any assert so it always lands) ---
    try
        fig = figure('Visible','off','Color','w','Position',[60 60 760 620]);
        ax  = axes('Parent',fig); set(ax,'Color','w'); hold(ax,'on');
        imagesc(ax, psf0); axis(ax,'image'); set(ax,'YDir','normal'); colormap(ax,hot);
        colorbar(ax); plot(ax, pc, pr, 'c+', 'MarkerSize',12, 'LineWidth',1.5);
        contour(ax, psf0, [0.5 0.5]*pk, 'c', 'LineWidth',1.0);
        title(ax, sprintf('PSF "%s": peak %.3g @(%d,%d), FWHM ax %.0f / lat %.0f um', ...
            fld, pk, pr, pc, fwhmAx_mm*1e3, fwhmLat_mm*1e3), 'Interpreter','none');
        xlabel(ax,'lateral (px)'); ylabel(ax,'axial (px)');
        local_check_disk_space(cfg.outputRoot, cfg.minFreeGB_write, 'psf validation png');
        exportgraphics(fig, fullfile(cfg.outputRoot,'psf_validation.png'), ...
            'Resolution',150, 'BackgroundColor','white');
        close(fig);
    catch ME
        fprintf('  WARN: psf_validation.png failed (%s)\n', ME.message);
    end

    % --- Hard validation (abort if insane) ---
    assert(all(isfinite(psf0(:))) && isreal(psf0), 'APR17:PSF', 'PSF not finite/real.');
    assert(pk > 0 && integ > 0, 'APR17:PSF', 'PSF not non-negative-dominant (peak/integral <= 0).');
    assert(pr > 0.10*Ny && pr < 0.90*Ny && pc > 0.10*Nx && pc < 0.90*Nx, ...
        'APR17:PSF', 'PSF peak not within the central region (row %d/%d, col %d/%d).', pr, Ny, pc, Nx);
    % single dominant peak: the above-half-max main lobe occupies a small
    % fraction of the image (a localized blob, not flat / noise / many peaks).
    fracMainLobe = nnz(psf0 > 0.5*pk) / numel(psf0);
    assert(fracMainLobe < 0.10, 'APR17:PSF', ...
        'PSF lacks a single dominant compact peak (above-half-max area = %.1f%% of image).', ...
        100*fracMainLobe);
    assert(fwhmAx_mm  > 0.01 && fwhmAx_mm  < 2.0, 'APR17:PSF', ...
        'PSF axial FWHM %.4f mm out of plausible range.', fwhmAx_mm);
    assert(fwhmLat_mm > 0.01 && fwhmLat_mm < 2.0, 'APR17:PSF', ...
        'PSF lateral FWHM %.4f mm out of plausible range.', fwhmLat_mm);

    % --- Compact odd-sided kernel around the peak (energy > 1e-3*max) ---
    m = psf0 > 1e-3*pk;
    [rr, ccx] = find(m);
    hr = max([pr-min(rr), max(rr)-pr]);
    hc = max([pc-min(ccx), max(ccx)-pc]);
    hr = min(hr, min(pr-1, Ny-pr));      % clip so the symmetric window stays in-bounds
    hc = min(hc, min(pc-1, Nx-pc));
    kernel = psf0(pr-hr:pr+hr, pc-hc:pc+hc);

    psf = struct('kernel', single(kernel), 'field', fld, ...
        'peak', pk, 'peakRow', pr, 'peakCol', pc, ...
        'kCenterRow', hr+1, 'kCenterCol', hc+1, ...
        'fwhmAx_mm', fwhmAx_mm, 'fwhmLat_mm', fwhmLat_mm, ...
        'fwhmAx_px', fwhmAx_px, 'fwhmLat_px', fwhmLat_px, ...
        'srcSize', [Ny Nx]);
    fprintf('  compact kernel: %dx%d (center at %d,%d)\n', ...
        size(kernel,1), size(kernel,2), psf.kCenterRow, psf.kCenterCol);
end

function w = local_fwhm_px(prof, pkIdx)
    prof = double(prof(:)); v = prof(pkIdx); half = 0.5*v;
    if v <= 0, w = NaN; return; end
    % walk left
    iL = pkIdx;
    while iL > 1 && prof(iL) > half, iL = iL - 1; end
    if prof(iL) > half, xL = iL;
    else, xL = iL + (half - prof(iL)) / max(prof(iL+1)-prof(iL), eps); end
    % walk right
    iR = pkIdx;
    while iR < numel(prof) && prof(iR) > half, iR = iR + 1; end
    if prof(iR) > half, xR = iR;
    else, xR = iR - (half - prof(iR)) / max(prof(iR-1)-prof(iR), eps); end
    w = xR - xL;
end

% Re-pad/re-center the compact kernel onto a [nZ x nX] grid, peak at the
% centered-PSF location (floor(nZ/2)+1, floor(nX/2)+1) for ifftshift.
function psfBlock = local_psf_to_grid(psf, g)
    psfBlock = zeros(g.nZ, g.nX, 'single');
    K = psf.kernel; [kh, kw] = size(K);
    cz = floor(g.nZ/2)+1; cx = floor(g.nX/2)+1;
    hr = (kh-1)/2; hc = (kw-1)/2;
    % crop kernel symmetrically if it exceeds the grid
    mr = min(hr, cz-1); mr = min(mr, g.nZ-cz);
    mc = min(hc, cx-1); mc = min(mc, g.nX-cx);
    kr = (psf.kCenterRow-mr):(psf.kCenterRow+mr);
    kc = (psf.kCenterCol-mc):(psf.kCenterCol+mc);
    psfBlock(cz-mr:cz+mr, cx-mc:cx+mc) = K(kr, kc);
end

% Derive a domain stack and SVD8-filter it once.
function IQf = local_domain_iqf(IQ_pos, IQ_neg, domain, cfg)
    switch domain
        case 'PI'
            if isempty(IQ_neg), D = IQ_pos; else, D = IQ_pos + IQ_neg; end
        case 'fundamental'
            if isempty(IQ_neg), D = IQ_pos; else, D = IQ_pos - IQ_neg; end
        case 'singlepol'
            D = IQ_pos;
        otherwise
            error('APR17:Domain', 'Unknown domain %s', domain);
    end
    IQf = svd_clutter_filter_rsvd(D, cfg.svdCutoff, [], cfg.useGPU, cfg.svdSeed);
end

% Per-block readout: amplitude (FISTA->refit->ROI sum) AND in-tube loc-rate
% (per-frame detect_microbubbles -> ROI mask count / nF; NO tracking) + density.
function r = local_block_readout(IQf, psfBlock, masks, g, thr, cfg)
    nROI = numel(masks);
    nF   = size(IQf, 3);

    % --- amplitude ---
    powerImg = mean(abs(IQf).^2, 3);
    [x_s, ~] = sushi_sparse_recovery(powerImg, psfBlock, cfg.lambda, 'fista', cfg.fista);
    [x_q, ~] = quasar_refit(powerImg, psfBlock, x_s, cfg.refit);
    x_q = gather(max(x_q, 0));
    x_s = gather(max(x_s, 0));
    ampQ = zeros(1, nROI); ampS = zeros(1, nROI);
    for rr = 1:nROI
        ampQ(rr) = sum(x_q(masks{rr}), 'all');
        ampS(rr) = sum(x_s(masks{rr}), 'all');
    end

    % --- in-tube localization RATE (no tracking) ---
    % NOTE: roiSize_px is unused by the 'fixed' detector (only minSep_mm drives NMS);
    % kept for parity with the ground-truth det struct.
    det = struct('method','fixed', 'fixedThresh',thr, 'roiSize_px',cfg.det.roiSize_px, ...
                 'minSep_mm', g.minSep_mm);
    allR = []; allC = [];
    for fr = 1:nF
        pix = detect_microbubbles(abs(IQf(:,:,fr)), det, g.dx, g.dz);   % [N x 2] = [row col]
        if isempty(pix), continue; end
        allR = [allR; pix(:,1)]; allC = [allC; pix(:,2)]; %#ok<AGROW>
    end
    locRate = zeros(1, nROI);
    if ~isempty(allR)
        lin = sub2ind([g.nZ g.nX], allR, allC);
        for rr = 1:nROI
            locRate(rr) = sum(masks{rr}(lin)) / nF;
        end
        dens = accumarray([allR allC], 1, [g.nZ g.nX]);
    else
        dens = zeros(g.nZ, g.nX);
    end

    r = struct('ampQ',ampQ, 'ampS',ampS, 'locRate',locRate, ...
               'locRate_fov', numel(allR)/nF, 'nF', nF, 'dens', dens);
end

% Per-domain detection threshold from the flowing-bg blocks.
%
% WHAT THIS ACTUALLY RETURNS: prctile(envPool, 99.9), the ceiling of its own candidate range.
% NOT a calibrated knee. The tolerance search below never succeeds and always falls through to
% the min(falseRate) branch, which selects the largest candidate, which is `hi`. Verified on
% all six (dataset x domain) combinations: every shipped threshold is within 2% of its own
% p99.9, and the measured Bg false-alarm rate at those thresholds is 12-18 loc/frame against
% the nominal tol of 0.02. See docs/FINDINGS_2026_07_14_threshold_ceiling_and_psf.md.
%
% WHY IT IS LEFT THIS WAY, DELIBERATELY:
%  (a) p99.9 is a defensible operating point, though not for the reason the code implies. It is
%      a bias/variance compromise: lower thr catches more dim bubbles but subtracts a bigger
%      pedestal (variance); higher thr subtracts less but is progressively amplitude-selected
%      (bias). The count is unbiased at ANY threshold because the noise floor is
%      concentration-independent, so the Bg pedestal subtraction is exact.
%  (b) REACHING the 0.02 tolerance would make things WORSE, not better. A false-alarm-free
%      threshold counts only the bright subset, and the bright fraction grows with
%      concentration (8.3% -> 51.1% across the jun23 ladder), fabricating a slope increase of
%      +0.364. The tolerance is the wrong design goal for this metric.
%  (c) The value is preserved exactly so all published results reproduce.
%
% The tolerance/knee machinery below is therefore ABANDONED-IN-PLACE, not repaired. Do not
% "fix" it by widening the candidate range. If you want a different operating point, change it
% deliberately and re-derive every downstream number, knowing (b).
function [thr, curve] = local_sweep_threshold(domain, bgCachePaths, roi, cfg)
    % Gather domain envelopes + per-block combinedTube masks + per-block grid for
    % each Bg block. Store the real envelope (abs) once (halves memory vs complex
    % IQf and avoids recomputing abs at every candidate threshold).
    iComb = find(strcmp(roi.names, 'combinedTube'), 1);
    envPool = [];
    bg = struct('env',{},'mask',{},'nF',{},'g',{});
    for b = 1:numel(bgCachePaths)
        cp = bgCachePaths{b};
        if ~exist(cp,'file'), continue; end
        L = load(cp, 'IQ_pos','IQ_neg','g');
        IQf = local_domain_iqf(L.IQ_pos, L.IQ_neg, domain, cfg);
        nF  = size(IQf, 3);
        if nF > cfg.thrSweep.nFrMax
            IQf = IQf(:,:,1:cfg.thrSweep.nFrMax); nF = cfg.thrSweep.nFrMax;
        end
        env  = abs(IQf); clear IQf;
        mask = local_poly2mask(roi.poly{iComb}, L.g.xGrid, L.g.zGrid);
        bg(end+1) = struct('env',env, 'mask',mask, 'nF',nF, 'g',L.g); %#ok<AGROW>
        envPool = [envPool; reshape(env(1:max(1,round(numel(env)/2e4)):end), [], 1)]; %#ok<AGROW>
    end
    if isempty(bg)
        warning('APR17:NoBg', 'No Bg caches for domain %s sweep; thr=NaN', domain);
        thr = NaN; curve = struct('thr',[],'falseRate',[]); return;
    end

    % candidate span from the Bg envelope percentiles (cfg.thrSweep.span='auto').
    lo  = prctile(envPool, 50);
    hi  = prctile(envPool, 99.9);
    cand = linspace(lo, hi, cfg.thrSweep.nThr);

    falseRate = zeros(1, cfg.thrSweep.nThr);
    for it = 1:cfg.thrSweep.nThr
        nIn = 0; nFr = 0;
        for b = 1:numel(bg)
            % SAME detector as the rungs (fixed; roiSize 5; minSep = lambda),
            % on THIS Bg block's own natural grid.
            det = struct('method','fixed', 'fixedThresh',cand(it), ...
                         'roiSize_px',cfg.det.roiSize_px, 'minSep_mm', bg(b).g.minSep_mm);
            for fr = 1:bg(b).nF
                pix = detect_microbubbles(bg(b).env(:,:,fr), det, bg(b).g.dx, bg(b).g.dz);
                if isempty(pix), continue; end
                lin = sub2ind(size(bg(b).mask), pix(:,1), pix(:,2));
                nIn = nIn + sum(bg(b).mask(lin));
            end
            nFr = nFr + bg(b).nF;
        end
        falseRate(it) = nIn / max(nFr, 1);
    end

    % Knee = lowest threshold ABOVE WHICH the clean condition PERSISTS (falseRate
    % <= tol for ALL higher candidates), not merely the first transient dip below
    % tol (Bg false-rate is a noisy estimate; a low-thr dip could be spurious).
    lastDirty = find(falseRate > cfg.thrSweep.tolLocPerFrame, 1, 'last');
    if isempty(lastDirty)
        okIdx = 1;                      % every candidate already clean
    else
        okIdx = lastDirty + 1;          % first candidate above the last dirty one
    end
    tolMet = true;
    if okIdx > numel(falseRate)
        % EXPECTED PATH, always taken. The candidate ceiling (prctile(envPool,99.9), whole-FOV)
        % sits below the in-tube peaks, so no candidate reaches tol and we return `hi` itself.
        % This is documented and deliberate, NOT a silent degradation: see the function header
        % and docs/FINDINGS_2026_07_14_threshold_ceiling_and_psf.md. The count stays unbiased
        % via the Bg pedestal subtraction (the noise floor is concentration-independent).
        tolMet = false;
        [~, okIdx] = min(falseRate);
        fprintf(['  [sweep %s] tol NOT met (expected): no candidate reaches %.3f loc/frame.\n' ...
                 '             thr = p99.9 CEILING of the candidate range, not a knee.\n' ...
                 '             Bg false rate at thr = %.2f loc/frame; the Bg PEDESTAL\n' ...
                 '             SUBTRACTION is what makes the count unbiased. Do NOT raise thr\n' ...
                 '             to chase tol: that amplitude-selects and inflates the slope.\n'], ...
            domain, cfg.thrSweep.tolLocPerFrame, falseRate(okIdx));
    end
    thr = cand(okIdx);
    % 'knee' is kept as a FIELD NAME for backward compatibility with saved readouts; it is the
    % p99.9 ceiling whenever tolMet is false, which is always. tolMet records which it is.
    curve = struct('thr',cand, 'falseRate',falseRate, 'knee',thr, 'tol',cfg.thrSweep.tolLocPerFrame, ...
                   'tolMet',tolMet, 'isCeiling',~tolMet);

    % --- ROC-like curve PNG (white) ---
    try
        fig = figure('Visible','off','Color','w','Position',[60 60 760 520]);
        ax = axes('Parent',fig); set(ax,'Color','w'); hold(ax,'on');
        semilogy(ax, cand, max(falseRate,1e-4), 'o-', 'Color',[0.1 0.3 0.7], 'LineWidth',1.6);
        yline(ax, cfg.thrSweep.tolLocPerFrame, '--', 'tol (never reached)', 'Color',[0.6 0 0]);
        if tolMet
            lbl = sprintf('knee=%.4g', thr); ttl = sprintf('Threshold sweep [%s]: knee=%.4g (tol=%.3f MET)', ...
                domain, thr, cfg.thrSweep.tolLocPerFrame);
        else
            lbl = sprintf('CEILING=%.4g (p99.9)', thr);
            ttl = sprintf(['Threshold [%s]: p99.9 CEILING=%.4g, NOT a knee (tol=%.3f unmet; Bg=%.1f loc/frame)\n' ...
                           'count stays unbiased via Bg pedestal subtraction; do NOT raise thr to chase tol'], ...
                domain, thr, cfg.thrSweep.tolLocPerFrame, falseRate(okIdx));
        end
        xline(ax, thr, '-', lbl, 'Color',[0 0.5 0]);
        xlabel(ax,'fixed envelope threshold'); ylabel(ax,'in-tube false loc/frame (Bg)');
        title(ax, ttl, 'Interpreter','none');
        local_check_disk_space(cfg.outputRoot, cfg.minFreeGB_write, ['thr sweep ' domain]);
        exportgraphics(fig, fullfile(cfg.outputRoot, sprintf('thr_sweep_%s.png',domain)), ...
            'Resolution',150, 'BackgroundColor','white');
        close(fig);
    catch ME
        fprintf('  WARN: thr_sweep_%s.png failed (%s)\n', domain, ME.message);
    end
    fprintf('  [sweep %s] knee thr = %.5g  (falseRate=%.4f /frame)\n', domain, thr, falseRate(okIdx));
end

% Orchestrate STEP 2: per-domain sweep -> readouts -> pedestals -> fits -> montages.
%   [VERBATIM sweep/readout/pedestal/aggregate core; profile edits: ds param,
%    maxBlk block dimension, label-based window names, first-bgFlow montage tag,
%    ds.name save names, ds passed to montage]
function [readout, betas] = local_step2_analyze(cfg, blocks, roi, psf, polCacheDir, ds)
    nDom = numel(cfg.domains);
    nB   = numel(blocks);
    nROI = numel(roi.names);
    rungMask = strcmp({blocks.kind},'rung');
    nRung  = max([blocks(rungMask).rungIdx]);
    maxBlk = max([blocks(rungMask).blk]);          % variable per dataset (apr17=3, jun23<=6)
    rungConc   = nan(1, nRung);
    rungLabels = cell(1, nRung);                    % for label-based window names
    for k = 1:nB
        if strcmp(blocks(k).kind,'rung')
            rungConc(blocks(k).rungIdx)   = blocks(k).conc;
            rungLabels{blocks(k).rungIdx} = blocks(k).rung;
        end
    end

    bgIdx   = find(strcmp({blocks.kind},'bgFlow') & ~[blocks.missing]);
    nofIdx  = find(strcmp({blocks.kind},'bgStatic') & ~[blocks.missing]);
    bgCachePaths = arrayfun(@(k) local_cache_path(polCacheDir, blocks(k).userLabel), bgIdx, 'uni', 0);
    firstBg = ds.bgFlow{1};                         % montage bubble-free tile label

    % Per-domain readout (cached; resume-safe within STEP 2).
    domData = cell(1, nDom);
    for d = 1:nDom
        dom = cfg.domains{d};
        cacheD = fullfile(cfg.outputRoot, sprintf('readout_%s.mat', dom));
        if exist(cacheD, 'file')
            fprintf('\n  ===== domain %s: cached readout, load =====\n', dom);
            Ld = load(cacheD, 'Rd'); domData{d} = Ld.Rd; continue;
        end
        fprintf('\n  ===== domain %s: sweep + per-block readout =====\n', dom);
        tD = tic;

        % 1) threshold (override iff finite, else sweep)
        ovr = cfg.thrFixed.(dom);
        if isfinite(ovr)
            thr = ovr; curve = struct('override',true,'thr',thr);
            fprintf('  [%s] using cfg.thrFixed override = %.5g\n', dom, thr);
        else
            [thr, curve] = local_sweep_threshold(dom, bgCachePaths, roi, cfg);
        end
        assert(isfinite(thr), 'APR17:Thr', 'domain %s threshold non-finite', dom);

        % 2) per-block readouts (load each cache ONCE)
        Rd = struct();
        Rd.domain   = dom;
        Rd.thr      = thr;
        Rd.curve    = curve;
        Rd.ampQ     = nan(nB, nROI);
        Rd.ampS     = nan(nB, nROI);
        Rd.locRate  = nan(nB, nROI);
        Rd.locFov   = nan(nB, 1);
        Rd.nF       = nan(nB, 1);
        Rd.fracInComb = nan(nB, 1);
        Rd.dens     = cell(nB, 1);     % stored only for montage blocks (b2 + first bgFlow)

        iComb = find(strcmp(roi.names,'combinedTube'),1);
        for k = 1:nB
            if blocks(k).missing, continue; end
            cp = local_cache_path(polCacheDir, blocks(k).userLabel);
            if ~exist(cp,'file')
                fprintf('  [%s] %s: no cache, skip\n', dom, blocks(k).userLabel); continue;
            end
            try
                L = load(cp, 'IQ_pos','IQ_neg','g');
                g = L.g;
                masks = cell(1, nROI);
                for rr = 1:nROI, masks{rr} = local_roi_mask(roi, rr, g); end
                psfBlock = local_psf_to_grid(psf, g);
                IQf = local_domain_iqf(L.IQ_pos, L.IQ_neg, dom, cfg);
                clear L;
                r = local_block_readout(IQf, psfBlock, masks, g, thr, cfg);
                clear IQf;
                Rd.ampQ(k,:)    = r.ampQ;
                Rd.ampS(k,:)    = r.ampS;
                Rd.locRate(k,:) = r.locRate;
                Rd.locFov(k)    = r.locRate_fov;
                Rd.nF(k)        = r.nF;
                Rd.fracInComb(k)= r.locRate(iComb) / max(r.locRate_fov, eps);
                isMont = (strcmp(blocks(k).kind,'rung') && blocks(k).blk==2) || ...
                         strcmp(blocks(k).userLabel, firstBg);
                if isMont, Rd.dens{k} = r.dens; end
                fprintf('  [%s] %-9s nF=%d locRate(comb)=%.3f ampQ(comb)=%.3g\n', ...
                    dom, blocks(k).userLabel, r.nF, r.locRate(iComb), r.ampQ(iComb));
            catch ME
                fprintf('  [%s] %s ERROR: %s\n', dom, blocks(k).userLabel, ME.message);
                local_print_stack(ME);
            end
        end

        % 3) pedestals (mean over the flowing-bg blocks)
        Rd.pedAmp = mean(Rd.ampQ(bgIdx,:), 1, 'omitnan');
        Rd.pedLoc = mean(Rd.locRate(bgIdx,:), 1, 'omitnan');

        try
            local_check_disk_space(cfg.outputRoot, cfg.minFreeGB_write, ['readout ' dom]);
            save(cacheD, 'Rd', '-v7.3');
        catch ME
            fprintf('  WARN: could not save %s (%s)\n', cacheD, ME.message);
        end
        fprintf('  [%s] done (%.1f min)\n', dom, toc(tD)/60);
        domData{d} = Rd;
    end

    % ---- Assemble tensors, aggregate, fit, montage ----
    readout = struct();
    readout.domains   = cfg.domains;
    readout.roiNames  = roi.names;
    readout.roi       = roi;                          % physical-polygon ROI provenance
    readout.rungConc  = rungConc;
    readout.rungLabels= rungLabels;
    readout.blocks    = blocks;
    readout.psf       = psf;
    readout.dens      = cell(1, nDom);                % per-block density (montage blocks only)
    readout.thrFixed  = nan(1, nDom);
    readout.ampQ      = nan(nDom, nRung, maxBlk, nROI);   % [dom rung blk roi]
    readout.ampS      = nan(nDom, nRung, maxBlk, nROI);
    readout.locRate   = nan(nDom, nRung, maxBlk, nROI);
    readout.nFrames   = nan(nRung, maxBlk);
    readout.pedAmp    = nan(nDom, nROI);
    readout.pedLoc    = nan(nDom, nROI);
    readout.noFlowAmp = nan(nDom, nROI);
    readout.noFlowLoc = nan(nDom, nROI);
    readout.bgAmp     = nan(nDom, numel(bgIdx), nROI);
    readout.bgLoc     = nan(nDom, numel(bgIdx), nROI);

    for d = 1:nDom
        Rd = domData{d};
        readout.thrFixed(d) = Rd.thr;
        readout.pedAmp(d,:) = Rd.pedAmp;
        readout.pedLoc(d,:) = Rd.pedLoc;
        readout.dens{d}     = Rd.dens;     % cell{nB}, montage blocks (b2 + first bgFlow) populated
        for k = 1:nB
            if blocks(k).missing, continue; end
            if strcmp(blocks(k).kind,'rung')
                iR = blocks(k).rungIdx; b = blocks(k).blk;
                readout.ampQ(d,iR,b,:)    = Rd.ampQ(k,:);
                readout.ampS(d,iR,b,:)    = Rd.ampS(k,:);
                readout.locRate(d,iR,b,:) = Rd.locRate(k,:);
                readout.nFrames(iR,b)     = Rd.nF(k);
            end
        end
        for j = 1:numel(bgIdx)
            readout.bgAmp(d,j,:) = Rd.ampQ(bgIdx(j),:);
            readout.bgLoc(d,j,:) = Rd.locRate(bgIdx(j),:);
        end
        if ~isempty(nofIdx)
            readout.noFlowAmp(d,:) = Rd.ampQ(nofIdx(1),:);
            readout.noFlowLoc(d,:) = Rd.locRate(nofIdx(1),:);
        end
    end

    % Per-rung aggregation (mean +/- SEM over the blocks).
    readout.ampQ_mean = squeeze(mean(readout.ampQ, 3, 'omitnan'));   % [dom rung roi]
    readout.ampS_mean = squeeze(mean(readout.ampS, 3, 'omitnan'));
    readout.loc_mean  = squeeze(mean(readout.locRate, 3, 'omitnan'));
    readout.ampQ_sem  = squeeze(local_sem(readout.ampQ, 3));
    readout.loc_sem   = squeeze(local_sem(readout.locRate, 3));

    % Bg-subtracted per-rung means (aggregate-then-subtract, ground-truth order).
    readout.ampQ_sub = nan(nDom, nRung, nROI);
    readout.loc_sub  = nan(nDom, nRung, nROI);
    for d = 1:nDom
        readout.ampQ_sub(d,:,:) = max(squeeze(readout.ampQ_mean(d,:,:)) - readout.pedAmp(d,:), 0);
        readout.loc_sub(d,:,:)  = max(squeeze(readout.loc_mean(d,:,:))  - readout.pedLoc(d,:), 0);
    end

    % ---- log-log beta fits per {domain x ROI x metric x window} ----
    betas = struct('rows', {{}}, 'lookup', struct());
    metrics = {'locRate','ampSUSHI','ampFISTA'};   % ampSUSHI=refit x_q (headline); ampFISTA=pre-refit x_s (diagnostic)
    for d = 1:nDom
        dom = cfg.domains{d};
        for rr = 1:nROI
            rn = roi.names{rr};
            for mi = 1:numel(metrics)
                met = metrics{mi};
                switch met
                    case 'locRate'
                        yraw = squeeze(readout.loc_mean(d,:,rr));
                        ysub = squeeze(readout.loc_sub(d,:,rr));
                    case 'ampSUSHI'    % refit output x_q = the headline "SUSHI amplitude"
                        yraw = squeeze(readout.ampQ_mean(d,:,rr));
                        ysub = squeeze(readout.ampQ_sub(d,:,rr));
                    case 'ampFISTA'    % pre-refit FISTA x_s = diagnostic ONLY
                        yraw = squeeze(readout.ampS_mean(d,:,rr));
                        ysub = yraw;   % reported RAW: NO Bg-pedestal subtraction (diagnostic)
                end
                for wi = 1:numel(cfg.fitWindows)
                    w = cfg.fitWindows{wi}; sel = w(1):w(2);
                    wname = sprintf('%s-%s', rungLabels{w(1)}, rungLabels{w(2)});  % label-based (C3-C7 / M1-M5)
                    fr = lf(rungConc(sel), ysub(sel), cfg.nBoot);
                    fwr = lf(rungConc(sel), yraw(sel), cfg.nBoot);
                    betas.rows{end+1} = struct('domain',dom,'roi',rn,'metric',met, ...
                        'window',wname, 'beta_sub',fr.beta,'ci_sub',fr.ci,'R2_sub',fr.R2, ...
                        'beta_raw',fwr.beta,'ci_raw',fwr.ci,'R2_raw',fwr.R2,'n',fr.n);
                    key = matlab.lang.makeValidName(sprintf('%s_%s_%s_%s', dom, rn, met, wname));
                    betas.lookup.(key) = betas.rows{end};
                end
            end
        end
    end

    % ---- print beta table ----
    fprintf('\n  ===== BETA TABLE (Bg-subtracted; raw in parentheses) =====\n');
    fprintf('  %-12s %-13s %-10s %-7s  %8s %16s %6s\n', ...
        'domain','roi','metric','window','beta','CI','R2');
    for i = 1:numel(betas.rows)
        R = betas.rows{i};
        fprintf('  %-12s %-13s %-10s %-7s  %8.3f [%.2f,%.2f] %5.2f  (raw %.3f)\n', ...
            R.domain, R.roi, R.metric, R.window, R.beta_sub, R.ci_sub(1), R.ci_sub(2), ...
            R.R2_sub, R.beta_raw);
    end

    % ---- montages (per domain) ----
    for d = 1:nDom
        try
            local_render_montage(cfg, domData{d}, blocks, roi, polCacheDir, ds);
        catch ME
            fprintf('  WARN: montage [%s] failed (%s)\n', cfg.domains{d}, ME.message);
            local_print_stack(ME);
        end
    end

    % ---- save (per-dataset names) ----
    try
        local_check_disk_space(cfg.outputRoot, cfg.minFreeGB_write, sprintf('readout_%s save', ds.name));
        save(fullfile(cfg.outputRoot, sprintf('readout_%s.mat', ds.name)), 'readout', '-v7.3');
    catch ME
        fprintf('  ERROR saving readout_%s.mat: %s\n', ds.name, ME.message); local_print_stack(ME);
    end
    try
        local_check_disk_space(cfg.outputRoot, cfg.minFreeGB_write, sprintf('betas_%s save', ds.name));
        save(fullfile(cfg.outputRoot, sprintf('betas_%s.mat', ds.name)), 'betas');
    catch ME
        fprintf('  ERROR saving betas_%s.mat: %s\n', ds.name, ME.message); local_print_stack(ME);
    end
end

% Localization-density montage (rep block/rung + first bgFlow) for one domain.
%   [generalized tile layout: tiles = blk==2 per rung + ds.bgFlow{1};
%    grid = ceil(sqrt(nRung+1)); titles from rungs.label/conc]
function local_render_montage(cfg, Rd, blocks, roi, polCacheDir, ds)
    dom = Rd.domain;
    iComb = find(strcmp(roi.names,'combinedTube'),1);
    iBg   = find(strcmp(roi.names,'background'),1);
    nRung  = numel(ds.rungs);
    nTiles = nRung + 1;
    gridN  = ceil(sqrt(nTiles));
    order = cell(1,nTiles);
    for iR = 1:nRung, order{iR} = sprintf(ds.blockFmt, ds.rungs(iR).label, 2); end  % blk==2 representative
    order{nTiles} = ds.bgFlow{1};
    tiles = cell(1,nTiles); labs = cell(1,nTiles);
    concL = nan(1,nTiles); fov = nan(1,nTiles); frac = nan(1,nTiles);
    for t = 1:nTiles
        k = find(strcmp({blocks.userLabel}, order{t}), 1);
        labs{t} = order{t};
        if isempty(k) || isempty(Rd.dens{k}), continue; end
        tiles{t} = Rd.dens{k};
        concL(t) = blocks(k).conc;
        fov(t)   = Rd.locFov(k);
        frac(t)  = Rd.fracInComb(k);
    end
    fig = figure('Visible','off','Color','w','Position',[30 30 1500 1050]);
    tl = tiledlayout(fig, gridN, gridN, 'Padding','compact','TileSpacing','compact');
    for t = 1:nTiles
        ax = nexttile(tl); set(ax,'Color','w');
        if isempty(tiles{t}), axis(ax,'off'); title(ax,[labs{t} ' (missing)']); continue; end
        k = find(strcmp({blocks.userLabel}, order{t}), 1);
        L = load(local_cache_path(polCacheDir, blocks(k).userLabel), 'g'); g = L.g;
        imagesc(ax, g.xGrid, g.zGrid, log10(tiles{t}+1)); colormap(ax,hot);
        axis(ax,'image'); set(ax,'YDir','normal','FontSize',8); hold(ax,'on');
        for rr = [iComb iBg]
            P = roi.poly{rr}; if isempty(P), continue; end
            cc = [0 1 1]; if rr==iBg, cc = [0 1 0]; end
            plot(ax, [P(:,1); P(1,1)], [P(:,2); P(1,2)], '-', 'Color', cc, 'LineWidth',1.2);
        end
        if isnan(concL(t))
            ttl = sprintf('%s (bubble-free) %.1f/fr in-tube %.0f%%', labs{t}, fov(t), 100*frac(t));
        else
            ttl = sprintf('%s %.2g MB/mL  %.1f/fr in-tube %.0f%%', labs{t}, concL(t), fov(t), 100*frac(t));
        end
        title(ax, ttl, 'FontSize',8.5, 'Interpreter','none');
        if t > gridN*(gridN-1), xlabel(ax,'x (mm)'); end
        if mod(t-1,gridN)==0,   ylabel(ax,'z (mm)'); end
    end
    title(tl, sprintf('Localization density [%s] (fixed thr=%.4g, log)  [cyan=combinedTube, green=background]', ...
        dom, Rd.thr), 'FontSize',11, 'Interpreter','none');
    outP = fullfile(cfg.outputRoot, sprintf('montage_localizations_%s.png', dom));
    local_check_disk_space(cfg.outputRoot, cfg.minFreeGB_write, ['montage ' dom]);
    exportgraphics(fig, outP, 'Resolution',150, 'BackgroundColor','white');
    close(fig);
    fprintf('  saved %s\n', outP);
end

% Log-log beta + 400-boot CI + R2 (verbatim ground-truth lf semantics).
function out = lf(x, y, nBoot)
    if nargin < 3 || isempty(nBoot), nBoot = 400; end
    x = x(:); y = y(:);
    v = isfinite(x) & isfinite(y) & x>0 & y>0; x = x(v); y = y(v);
    out.n = numel(x);
    if numel(x) < 3
        out.beta = NaN; out.ci = [NaN NaN]; out.R2 = NaN; return;
    end
    X = log10(x); Y = log10(y); p = polyfit(X, Y, 1); yh = polyval(p, X);
    out.R2 = 1 - sum((Y-yh).^2) / max(sum((Y-mean(Y)).^2), eps);
    n = numel(X); bs = nan(nBoot,1);
    for b = 1:nBoot
        idx = randi(n, n, 1); pb = polyfit(X(idx), Y(idx), 1); bs(b) = pb(1);
    end
    out.beta = p(1); out.ci = quantile(bs, [.025 .975]);
end

% Per-dataset validation/report. apr17: SOFT PASS/CHECK vs ds.target. jun23
% (ds.target empty): context betas only. Always prints combinedTube/tubeL/tubeR/
% full over the headline window + ds.notes. Provenance of the apr17 1.19/0.94
% targets: the 06-15 LEFT/loose single-tube ROI at a HARDCODED thr=480; this
% pipeline re-sweeps + uses DRAWN ROIs, so a near-miss = ROI/threshold drift,
% NOT a recipe regression -- it REPORTS (PASS/CHECK), never errors.
function local_report_validation(betas, ds)
    winName = sprintf('%s-%s', ds.rungs(ds.headlineWin(1)).label, ds.rungs(ds.headlineWin(2)).label);
    fprintf('\n  ===== VALIDATION [%s] (PI, %s, Bg-subtracted) =====\n', ds.name, winName);
    if ~isempty(ds.target)
        t = ds.target;
        fprintf(['  SOFT target provenance: 1.19/0.94 are from the 06-15 LEFT/loose single-tube\n' ...
                 '  ROI at hardcoded thr=480; this run re-sweeps + uses drawn ROIs, so a near-miss\n' ...
                 '  is ROI/threshold drift (+/- %.2f), never a hard fail. tubeL ~ the GT loose box.\n'], t.tol);
        local_report_beta(betas, 'PI', t.roi, 'locRate',  t.window, t.locBeta, t.tol);
        local_report_beta(betas, 'PI', t.roi, 'ampSUSHI', t.window, t.ampBeta, t.tol);
    else
        fprintf('  No locked target for %s; consistency report only (see dual-tube report).\n', ds.name);
    end

    fprintf('\n  ----- context (PI, %s, no gate) -----\n', winName);
    for rn = {'combinedTube','tubeL','tubeR','full'}
        local_report_beta(betas, 'PI', rn{1}, 'locRate',  winName, NaN, NaN);
        local_report_beta(betas, 'PI', rn{1}, 'ampSUSHI', winName, NaN, NaN);
    end
    if isfield(ds,'notes') && ~isempty(ds.notes)
        fprintf('\n  [note:%s] %s\n', ds.name, ds.notes);
    end
end

% Print one beta row; if target is finite, append a soft PASS/CHECK verdict.
function local_report_beta(betas, dom, roiName, met, window, target, tol)
    key = matlab.lang.makeValidName(sprintf('%s_%s_%s_%s', dom, roiName, met, window));
    if ~isfield(betas.lookup, key)
        fprintf('  %-3s %-13s %-9s %-6s : KEY MISSING\n', dom, roiName, met, window);
        return;
    end
    R = betas.lookup.(key);
    if isfinite(target)
        verdict = ternary(abs(R.beta_sub - target) <= tol, 'PASS', 'CHECK');
        fprintf('  %-3s %-13s %-9s %-6s : beta=%.3f CI[%.2f,%.2f] (target %.2f +/-%.2f) -> %s\n', ...
            dom, roiName, met, window, R.beta_sub, R.ci_sub(1), R.ci_sub(2), target, tol, verdict);
    else
        fprintf('  %-3s %-13s %-9s %-6s : beta=%.3f CI[%.2f,%.2f] R2=%.2f (raw %.3f)\n', ...
            dom, roiName, met, window, R.beta_sub, R.ci_sub(1), R.ci_sub(2), R.R2_sub, R.beta_raw);
    end
end

% Fetch a beta row struct from betas.lookup (or a NaN stub if absent).
function R = local_lookup_beta(betas, dom, roiName, met, window)
    key = matlab.lang.makeValidName(sprintf('%s_%s_%s_%s', dom, roiName, met, window));
    if isfield(betas.lookup, key)
        R = betas.lookup.(key);
    else
        R = struct('beta_sub',NaN,'ci_sub',[NaN NaN],'R2_sub',NaN,'beta_raw',NaN);
    end
end

% CI-overlap verdict between two beta rows.
function v = local_ci_verdict(a, b)
    if ~all(isfinite([a.beta_sub, b.beta_sub, a.ci_sub(:)', b.ci_sub(:)']))
        v = 'N/A (insufficient rungs)'; return;
    end
    if ~(a.ci_sub(2) < b.ci_sub(1) || b.ci_sub(2) < a.ci_sub(1))
        v = 'PASS (CIs overlap)';
    else
        v = 'CHECK (CIs disjoint)';
    end
end

% ---------------------------------------------------------------------
%   Dual-tube consistency (NEW): same per-rung concentration in BOTH tubes
%   => same beta. Quantifies tubeL-vs-tubeR agreement from the already-
%   computed readout/betas (no re-processing). Profile-driven expected ratio
%   (apr17 >1 size-driven; jun23 ~1). Writes <ds>/dualtube_consistency.{txt,png}.
% ---------------------------------------------------------------------
function local_dualtube_consistency(readout, betas, ds)
    minFreeGB = 3;   % == cfg.minFreeGB_write (cfg not threaded per spec signature)
    if ds.nTubes < 2
        fprintf('  [dualtube:%s] nTubes < 2; skipping.\n', ds.name); return;
    end
    winName = sprintf('%s-%s', ds.rungs(ds.headlineWin(1)).label, ds.rungs(ds.headlineWin(2)).label);
    dPI = find(strcmp(readout.domains,'PI'),1);
    iL  = find(strcmp(readout.roiNames,'tubeL'),1);
    iR  = find(strcmp(readout.roiNames,'tubeR'),1);
    if isempty(dPI) || isempty(iL) || isempty(iR)
        fprintf('  [dualtube:%s] PI/tubeL/tubeR not all present; skipping.\n', ds.name); return;
    end

    % ---- betas (PI, tubeL vs tubeR) over the headline window ----
    bL_loc = local_lookup_beta(betas,'PI','tubeL','locRate', winName);
    bR_loc = local_lookup_beta(betas,'PI','tubeR','locRate', winName);
    bL_amp = local_lookup_beta(betas,'PI','tubeL','ampSUSHI',winName);
    bR_amp = local_lookup_beta(betas,'PI','tubeR','ampSUSHI',winName);
    vLoc = local_ci_verdict(bL_loc, bR_loc);
    vAmp = local_ci_verdict(bL_amp, bR_amp);

    % ---- per-rung tubeL/tubeR ratios ----
    conc = readout.rungConc;
    locL = squeeze(readout.loc_sub(dPI,:,iL));  locR = squeeze(readout.loc_sub(dPI,:,iR));   % Bg-sub (match betas)
    ampL = squeeze(readout.ampQ_sub(dPI,:,iL)); ampR = squeeze(readout.ampQ_sub(dPI,:,iR));
    locLm = squeeze(readout.loc_mean(dPI,:,iL)); locRm = squeeze(readout.loc_mean(dPI,:,iR)); % raw means for stable ratio
    ampLm = squeeze(readout.ampQ_mean(dPI,:,iL)); ampRm = squeeze(readout.ampQ_mean(dPI,:,iR));
    ratLoc = locLm ./ locRm;  ratAmp = ampLm ./ ampRm;

    withdraw = isfield(ds.tubeNotes,'R') && ~isempty(regexpi(ds.tubeNotes.R, 'withdraw', 'once'));
    deviated = startsWith(vLoc,'CHECK') || startsWith(vAmp,'CHECK');

    % ---- text report ----
    txt = fullfile(ds.outDir, 'dualtube_consistency.txt');
    fid = fopen(txt, 'w');
    if fid > 0
        cu = onCleanup(@() fclose(fid));
        fprintf(fid, 'DUAL-TUBE CONSISTENCY [%s]   PI domain   window %s\n', ds.name, winName);
        fprintf(fid, 'expected tubeL/tubeR ratio: %s\n', ds.tubeNotes.expectRatio);
        fprintf(fid, 'LEFT  tube: %s\n', ds.tubeNotes.L);
        fprintf(fid, 'RIGHT tube: %s\n\n', ds.tubeNotes.R);
        fprintf(fid, '%-9s %20s %20s   %s\n', 'metric','beta_L [CI]','beta_R [CI]','CI-overlap');
        fprintf(fid, '%-9s %7.3f [%.2f,%.2f] %7.3f [%.2f,%.2f]   %s\n', 'locRate', ...
            bL_loc.beta_sub, bL_loc.ci_sub(1), bL_loc.ci_sub(2), ...
            bR_loc.beta_sub, bR_loc.ci_sub(1), bR_loc.ci_sub(2), vLoc);
        fprintf(fid, '%-9s %7.3f [%.2f,%.2f] %7.3f [%.2f,%.2f]   %s\n', 'ampSUSHI', ...
            bL_amp.beta_sub, bL_amp.ci_sub(1), bL_amp.ci_sub(2), ...
            bR_amp.beta_sub, bR_amp.ci_sub(1), bR_amp.ci_sub(2), vAmp);
        fprintf(fid, '\nper-rung tubeL/tubeR ratio (PI, raw aggregated means):\n');
        fprintf(fid, '%-6s %12s %10s %10s\n', 'rung','conc(MB/mL)','loc_ratio','amp_ratio');
        for iRu = 1:numel(ds.rungs)
            fprintf(fid, '%-6s %12.3g %10.3f %10.3f\n', ds.rungs(iRu).label, conc(iRu), ratLoc(iRu), ratAmp(iRu));
        end
        if withdraw && deviated
            fprintf(fid, ['\nNOTE: right-tube deviation is consistent with the logged WITHDRAW ' ...
                'episode;\n      use LEFT (clean reference) for the headline. No block auto-excluded.\n']);
        elseif withdraw
            fprintf(fid, ['\nNOTE: RIGHT tube had a brief logged withdraw episode but betas agree;\n' ...
                '      LEFT remains the reference. No block auto-excluded.\n']);
        end
        clear cu;
    else
        fprintf('  [dualtube:%s] WARN: could not open %s for writing.\n', ds.name, txt);
    end
    fprintf('  [dualtube:%s] loc %s ; amp %s -> %s\n', ds.name, vLoc, vAmp, txt);

    % ---- figure: overlaid tubeL/tubeR log-log + ratio-vs-rung (white) ----
    try
        fig = figure('Visible','off','Color','w','Position',[40 40 1500 460]);
        tl = tiledlayout(fig,1,3,'Padding','compact','TileSpacing','compact');

        ax1 = nexttile(tl); set(ax1,'Color','w'); hold(ax1,'on');
        loglog(ax1, conc, max(locL,eps), 'o-', 'Color',[0.10 0.45 0.80], 'LineWidth',1.6, 'DisplayName','tubeL');
        loglog(ax1, conc, max(locR,eps), 's--','Color',[0.85 0.33 0.10], 'LineWidth',1.6, 'DisplayName','tubeR');
        set(ax1,'XScale','log','YScale','log'); grid(ax1,'on');
        xlabel(ax1,'concentration (MB/mL)'); ylabel(ax1,'loc-rate (Bg-sub)');
        title(ax1, sprintf('loc: bL=%.2f bR=%.2f', bL_loc.beta_sub, bR_loc.beta_sub));
        legend(ax1,'Location','best');

        ax2 = nexttile(tl); set(ax2,'Color','w'); hold(ax2,'on');
        loglog(ax2, conc, max(ampL,eps), 'o-', 'Color',[0.10 0.45 0.80], 'LineWidth',1.6, 'DisplayName','tubeL');
        loglog(ax2, conc, max(ampR,eps), 's--','Color',[0.85 0.33 0.10], 'LineWidth',1.6, 'DisplayName','tubeR');
        set(ax2,'XScale','log','YScale','log'); grid(ax2,'on');
        xlabel(ax2,'concentration (MB/mL)'); ylabel(ax2,'SUSHI amp (Bg-sub)');
        title(ax2, sprintf('amp: bL=%.2f bR=%.2f', bL_amp.beta_sub, bR_amp.beta_sub));
        legend(ax2,'Location','best');

        ax3 = nexttile(tl); set(ax3,'Color','w'); hold(ax3,'on');
        plot(ax3, 1:numel(ratLoc), ratLoc, 'o-', 'Color',[0.10 0.45 0.80],'LineWidth',1.6,'DisplayName','loc L/R');
        plot(ax3, 1:numel(ratAmp), ratAmp, 's--','Color',[0.85 0.33 0.10],'LineWidth',1.6,'DisplayName','amp L/R');
        yline(ax3, 1, ':', 'Color',[0.4 0.4 0.4]);
        sel = ds.headlineWin(1):ds.headlineWin(2);
        xline(ax3, sel(1),  '-', 'Color',[0 0.5 0]);
        xline(ax3, sel(end),'-', 'Color',[0 0.5 0]);
        xlabel(ax3,'rung index'); ylabel(ax3,'tubeL / tubeR ratio'); grid(ax3,'on');
        title(ax3, sprintf('ratio (expect %s)', ds.tubeNotes.expectRatio), 'Interpreter','none');
        legend(ax3,'Location','best');

        title(tl, sprintf('Dual-tube consistency [%s]  (PI, window %s)', ds.name, winName), 'Interpreter','none');
        outP = fullfile(ds.outDir, 'dualtube_consistency.png');
        local_check_disk_space(ds.outDir, minFreeGB, ['dualtube ' ds.name]);
        exportgraphics(fig, outP, 'Resolution',150, 'BackgroundColor','white');
        close(fig);
        fprintf('  [dualtube:%s] saved %s\n', ds.name, outP);
    catch ME
        fprintf('  [dualtube:%s] figure failed (%s)\n', ds.name, ME.message);
    end
end

function s = ternary(cond, a, b)
    if cond, s = a; else, s = b; end
end

function se = local_sem(M, dim)
    n  = sum(isfinite(M), dim);
    sd = std(M, 0, dim, 'omitnan');
    se = sd ./ sqrt(max(n, 1));
end

% =====================================================================
%   REUSED hardening + skeleton helpers (from BATCH_REANALYSIS_2DATASET.m)
% =====================================================================
function seriesDir = local_find_series(studyRoot)
    if exist(fullfile(studyRoot, 'Series 1'), 'dir')
        seriesDir = fullfile(studyRoot, 'Series 1'); return;
    end
    subdirs = dir(studyRoot);
    for k = 1:numel(subdirs)
        if subdirs(k).isdir && ~startsWith(subdirs(k).name, '.')
            candidate = fullfile(studyRoot, subdirs(k).name, 'Series 1');
            if exist(candidate, 'dir'), seriesDir = candidate; return; end
        end
    end
    if ~isempty(dir(fullfile(studyRoot, ['*' '.vada.xml'])))
        seriesDir = studyRoot; return;
    end
    error('Cannot find Series 1 in %s', studyRoot);
end

function resolvedName = local_resolve_blockname(seriesDir, userBlockName)
% Word-boundary-safe label -> VADA basename (C5b1 != C5NoFlow; C5NoFlow exact).
    if ~isempty(dir(fullfile(seriesDir, [userBlockName '.vada.xml'])))
        resolvedName = userBlockName; return;
    end
    xmlFiles = dir(fullfile(seriesDir, '*.vada.xml'));
    matches = {};
    for k = 1:numel(xmlFiles)
        basename = regexprep(xmlFiles(k).name, '\.vada\.xml$', '');
        if contains(basename, userBlockName, 'IgnoreCase', true)
            matches{end+1} = basename; %#ok<AGROW>
        end
    end
    pat = ['(^|[^a-zA-Z0-9])' regexptranslate('escape', userBlockName) '($|[^a-zA-Z0-9])'];
    if isscalar(matches)
        if ~isempty(regexp(matches{1}, pat, 'once', 'ignorecase'))
            resolvedName = matches{1};
        else
            error('Block "%s" not found in %s', userBlockName, seriesDir);
        end
    elseif numel(matches) > 1
        bestIdx = 1;
        for k = 1:numel(matches)
            if ~isempty(regexp(matches{k}, pat, 'once', 'ignorecase')), bestIdx = k; break; end
        end
        resolvedName = matches{bestIdx};
    else
        error('Block "%s" not found in %s', userBlockName, seriesDir);
    end
end

function t = accum_time(t, fieldName, dt)
    if isfield(t, fieldName), t.(fieldName) = t.(fieldName) + dt; else, t.(fieldName) = dt; end
end

function freeGB = local_free_space_gb(pathOnDrive)
    freeGB = NaN;
    try
        jFile = java.io.File(pathOnDrive);
        bytes = jFile.getFreeSpace();
        if bytes > 0, freeGB = double(bytes) / 1e9; end
    catch
    end
end

function local_check_disk_space(pathOnDrive, minFreeGB, label)
    freeGB = local_free_space_gb(pathOnDrive);
    if isnan(freeGB)
        fprintf('  [disk] WARNING: could not query free space for %s (%s)\n', pathOnDrive, label);
        return;
    end
    fprintf('  [disk] %s: %.1f GB free on %s\n', label, freeGB, pathOnDrive);
    if freeGB < minFreeGB
        error('APR17:DiskLow', ['DISK SPACE GUARD: only %.1f GB free on drive containing %s ' ...
            '(need >= %.1f GB for %s). Aborting BEFORE write.'], freeGB, pathOnDrive, minFreeGB, label);
    end
end

function local_print_stack(ME)
    for es = 1:numel(ME.stack)
        fprintf('    at %s (line %d)\n', ME.stack(es).name, ME.stack(es).line);
    end
end
