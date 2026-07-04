%% BATCH_RUN.m
% Master batch orchestration: pipeline processing + presentation outputs.
%
% For each dataset:
%   1. Run LAT-ULM pipeline
%   2. Run QUASAR pipeline
%   3. Generate presentation figures (three-panel, FWHM, velocity profile)
%   4. Generate animation
%   5. Generate LAT-ULM vs QUASAR comparison
%
% USAGE:
%   - Configure datasets in the DATASETS section below
%   - Run the script; it will list all jobs and prompt for selection
%   - Overnight batch: select 'all' and let it run
%
% Author: Eli Wirth-Apley / Sun Lab, Northeastern University, March 2026

clearvars; close all; clc;

%% ========================================================================
%  PATHS
%  ========================================================================
baseDir         = 'C:\path\to\ULM3';
vadaScriptsPath = 'C:\path\to\VevoF2\RFExampleScriptsF2\Matlab';
addpath(genpath(vadaScriptsPath));
pipelineDir = fileparts(mfilename('fullpath'));
if isempty(pipelineDir), pipelineDir = pwd; end
addpath(pipelineDir);

%% ========================================================================
%  DATASETS
%  ========================================================================
% Each dataset defines: name, data folder, block filenames, background,
% expected velocities (for presentation figures), and which pipelines to run.
%
% IMPORTANT: Edit block lists below to match your actual filenames.
%            Use VADA_organize.m or dir() to discover block names.

DS = struct();
nDS = 0;

% --- Dataset 1: Phantom, initial speed ---
nDS = nDS+1;
DS(nDS).name        = 'Phantom Initial Speed';
DS(nDS).shortName   = 'phantom_initial';
DS(nDS).dataFolder  = fullfile(baseDir, 'Primary_Low');
DS(nDS).bgFile      = '2026-03-10-10-56-09_Waterbath_2_bg-2026-02-25-22-13-05';          % Background filename (no extension). '' = none
DS(nDS).bgFolder    = fullfile(baseDir, 'Primary_Low');
DS(nDS).blocks      = {'2026-03-10-10-56-09_2026-02-25-22-16-10', '2026-03-10-10-56-09_2026-02-25-22-18-22', '2026-03-10-10-56-09_2026-02-25-22-19-50', '2026-03-10-10-56-09_2026-02-25-22-21-30', '2026-03-10-10-56-09_2026-02-25-22-22-56', '2026-03-10-10-56-09_2026-02-25-22-24-53', '2026-03-10-10-56-09_2026-02-25-22-26-24', '2026-03-10-10-56-09_2026-02-25-22-28-06', '2026-03-10-10-56-09_2026-02-25-22-29-54', '2026-03-10-10-56-09_2026-02-25-22-31-28'};          % FILL IN: block filenames (no .vada extension)
DS(nDS).quasarBlocks = {'2026-03-10-10-56-09_2026-02-25-22-19-50', '2026-03-10-10-56-09_2026-02-25-22-21-30'};         % FILL IN: 1-3 blocks for QUASAR (subset of blocks above)
DS(nDS).outputBase  = fullfile(baseDir, 'Primary_Low', 'Results_InitialSpeed');
DS(nDS).runLATULM   = true;
DS(nDS).runQUASAR   = true;
DS(nDS).runPresentation = true;
DS(nDS).runAnimation    = true;
DS(nDS).runComparison   = true;

% Tube info for presentation figures
DS(nDS).tubes = struct( ...
    'name',       {'Small tube', 'Large tube'}, ...
    'ID_mm',      {0.31,          0.51}, ...
    'rate_mL_min',{0.08,          0.16}, ...
    'approxLateral_mm', {-.6,    .5});

% --- Dataset 2: Phantom, doubled speed ---
nDS = nDS+1;
DS(nDS).name        = 'Phantom Doubled Speed';
DS(nDS).shortName   = 'phantom_doubled';
DS(nDS).dataFolder  = fullfile(baseDir, 'Primary_Low');
DS(nDS).bgFile      = '2026-03-10-10-56-09_Waterbath_2_bg-2026-02-25-22-13-05';
DS(nDS).bgFolder    = fullfile(baseDir, 'Primary_Low');
DS(nDS).blocks      = {'2026-03-10-10-56-09_2026-02-25-22-32-58', '2026-03-10-10-56-09_2026-02-25-22-34-28', '2026-03-10-10-56-09_2026-02-25-22-36-13', '2026-03-10-10-56-09_2026-02-25-22-37-48'};          % FILL IN: last few block filenames (doubled pump rate)
DS(nDS).quasarBlocks = {'2026-03-10-10-56-09_2026-02-25-22-34-28', '2026-03-10-10-56-09_2026-02-25-22-36-13'};         % FILL IN: 1-3 blocks for QUASAR
DS(nDS).outputBase  = fullfile(baseDir, 'Primary_Low', 'Results_DoubledSpeed');
DS(nDS).runLATULM   = true;
DS(nDS).runQUASAR   = true;
DS(nDS).runPresentation = true;
DS(nDS).runAnimation    = true;
DS(nDS).runComparison   = true;

DS(nDS).tubes = struct( ...
    'name',       {'Small tube', 'Large tube'}, ...
    'ID_mm',      {0.31,          0.51}, ...
    'rate_mL_min',{0.16,          0.32}, ...  % Doubled
    'approxLateral_mm', {-.6,    .5});

% --- Dataset 3: Macrophage + Microbubble, initial speed ---
nDS = nDS+1;
DS(nDS).name        = 'Mac+MB Initial Speed';
DS(nDS).shortName   = 'mac_mb_initial';
DS(nDS).dataFolder  = 'C:\path\to\data\Series1';          % FILL IN: path to mac+MB VADA data
DS(nDS).bgFile      = '2026-03-11-15-26-17_Final bg vada-2026-03-11-14-24-23';          % FILL IN
DS(nDS).bgFolder    = 'C:\path\to\data\Series1';          % FILL IN
DS(nDS).blocks      = {'block_1', 'block_2', 'block_3', 'block_4', 'block_5', 'block_6', 'block_7', 'block_8', 'block_9', 'block_10'};          % FILL IN: all blocks EXCEPT the last 4
DS(nDS).quasarBlocks = {'block_2', 'block_3'};         % FILL IN: 1-3 blocks for QUASAR
DS(nDS).outputBase  = 'C:\path\to\data\Series1\Results_Mac&MB_Initial_Speed';          % FILL IN
DS(nDS).runLATULM   = true;
DS(nDS).runQUASAR   = true;
DS(nDS).runPresentation = true;
DS(nDS).runAnimation    = true;
DS(nDS).runComparison   = true;

DS(nDS).tubes = struct( ...        % Updated from mac+MB three_panel.png
    'name',       {'Small tube', 'Large tube'}, ...
    'ID_mm',      {0.31,          0.51}, ...
    'rate_mL_min',{0.08,          0.16}, ...
    'approxLateral_mm', {-0.25,    0.15});

% Mac+MB overrides: noisier echoes from macrophage-loaded MBs need
% moderately tighter clutter rejection and tracking filters
DS(nDS).configOverride.svd.cutoffLow        = 2;    % same as common — waterbath phantom, no tissue clutter
DS(nDS).configOverride.det.threshold        = 7;    % common 6
DS(nDS).configOverride.track.minTrackLength = 25;   % common 20
DS(nDS).configOverride.track.maxGapFrames   = 5;    % common 5 (same)
DS(nDS).configOverride.track.maxDisp_mm     = 0.050;% common 0.07
DS(nDS).configOverride.postproc.minSpeed_mm_s = 2.0;% common 1.5

% --- Dataset 4: Macrophage + Microbubble, doubled speed (last 4 blocks) ---
nDS = nDS+1;
DS(nDS).name        = 'Mac+MB Doubled Speed';
DS(nDS).shortName   = 'mac_mb_doubled';
DS(nDS).dataFolder  = 'C:\path\to\data\Series1';          % FILL IN: path to mac+MB VADA data
DS(nDS).bgFile      = '2026-03-11-15-26-17_Final bg vada-2026-03-11-14-24-23';          % FILL IN
DS(nDS).bgFolder    = 'C:\path\to\data\Series1';          % FILL IN
DS(nDS).blocks      = {'block_11', 'block_12', 'block_13', 'block_14'};          % FILL IN: last 4 block filenames only
DS(nDS).quasarBlocks = {'block_13', 'block_14'};         % FILL IN: 1-2 blocks for QUASAR
DS(nDS).outputBase  = 'C:\path\to\data\Series1\Results_Mac&MB_Double_Speed';          % FILL IN
DS(nDS).runLATULM   = true;
DS(nDS).runQUASAR   = true;
DS(nDS).runPresentation = true;
DS(nDS).runAnimation    = true;
DS(nDS).runComparison   = true;

DS(nDS).tubes = struct( ...        % Same geometry as DS(3)
    'name',       {'Small tube', 'Large tube'}, ...
    'ID_mm',      {0.31,          0.51}, ...
    'rate_mL_min',{0.16,          0.32}, ...  % Doubled
    'approxLateral_mm', {-0.25,    0.15});

% Same Mac+MB overrides as DS(3)
DS(nDS).configOverride.svd.cutoffLow        = 2;
DS(nDS).configOverride.det.threshold        = 7;
DS(nDS).configOverride.track.minTrackLength = 25;
DS(nDS).configOverride.track.maxGapFrames   = 5;
DS(nDS).configOverride.track.maxDisp_mm     = 0.050;
DS(nDS).configOverride.postproc.minSpeed_mm_s = 2.0;

%% ========================================================================
%  COMMON CONFIG (tuned for parallel flow phantom)
%  ========================================================================
C.useGPU = true;
C.numAngles = 3; C.numPolarities = 2; C.eventsPerFrame = 6;
C.zeroOnly = true; C.blankSteering = true; C.blankMargin = 1.5;
C.blankVoltage = []; C.txVoltage = [];
C.sosOverride = []; C.pitchOverride = 0.300; C.modeName = '.vada';
C.vadaScriptsPath = vadaScriptsPath; C.chunkSize = 400;

% Beamforming grid
C.bf.xRange = [-5, 5]; C.bf.zRange = [12, 24];
C.bf.dx = 0.025; C.bf.dz = 0.010;

% SVD
C.svd.cutoffLow = 2; C.svd.cutoffHigh = [];

% Detection (moderate for 5e5 — between old tight=7 and loose=5)
C.det.method = 'threshold'; C.det.threshold = 6;
C.det.minSep_mm = 0.100; C.det.roiSize_px = 7;

% Tracking (moderate — reject short noise tracks but accept 5e5 flow)
C.track.maxDisp_mm = 0.070;
C.track.maxGapFrames = 5;
C.track.minTrackLength = 20;
C.track.kalman.processNoise = 0.0005;
C.track.kalman.measNoise = 0.080;

% Super-resolution
C.sr.pixelSize_um = 5;

% SUSHI/QUASAR
C.ensemble.size = 150; C.ensemble.overlap = 0;
C.sushi.srFactor = 5; C.sushi.lambda = 0.1;
C.sushi.maxIter = 100; C.sushi.method = 'fista'; C.sushi.nonNeg = true;
C.quasar.enable = true; C.quasar.supportThresh = 0; C.quasar.maxIterCG = 50;
C.velocity.enable = true;

% Post-processing
C.postproc.minSpeed_mm_s = 1.5;     % Filter stationary/wall-stuck tracks
C.postproc.srPixel_um    = 25;      % Presentation figure SR pixel
C.postproc.velMax        = 30;      % Velocity colorbar max [mm/s]
C.postproc.anim_duration_s = 60;    % Animation clip length
C.postproc.anim_fps     = 30;
C.postproc.anim_smoothWindow = 5;
C.postproc.fwhm_depths  = [16.0, 16.5, 17.0, 17.5, 18.0];
C.postproc.velProfile_depth_mm = 17.0;

%% ========================================================================
%  BUILD JOB LIST
%  ========================================================================
fprintf('================================================================\n');
fprintf('  BATCH ORCHESTRATOR\n');
fprintf('================================================================\n\n');

jobs = struct('dsIdx',{}, 'pipeline',{}, 'name',{}, 'outDir',{}, 'blocks',{});

for d = 1:nDS
    ds = DS(d);
    if isempty(ds.blocks)
        fprintf('  %-25s -> SKIPPED (no blocks defined)\n', ds.name);
        continue;
    end
    if isempty(ds.dataFolder) || ~exist(ds.dataFolder, 'dir')
        fprintf('  %-25s -> SKIPPED (data folder missing)\n', ds.name);
        continue;
    end
    
    if ds.runLATULM
        jobs(end+1).dsIdx = d;
        jobs(end).pipeline = 'latulm';
        jobs(end).name = sprintf('%s — LAT-ULM', ds.name);
        jobs(end).outDir = fullfile(ds.outputBase, 'LAT_ULM_Results');
        jobs(end).blocks = ds.blocks;  % All blocks for LAT-ULM
    end
    if ds.runQUASAR
        % Use quasarBlocks if specified, otherwise default to first 2 blocks
        if ~isempty(ds.quasarBlocks)
            qBlocks = ds.quasarBlocks;
        else
            qBlocks = ds.blocks(1:min(2, numel(ds.blocks)));
            fprintf('  NOTE: %s QUASAR — no quasarBlocks set, defaulting to first %d block(s)\n', ...
                ds.name, numel(qBlocks));
        end
        jobs(end+1).dsIdx = d;
        jobs(end).pipeline = 'quasar';
        jobs(end).name = sprintf('%s — QUASAR', ds.name);
        jobs(end).outDir = fullfile(ds.outputBase, 'QUASAR_Results');
        jobs(end).blocks = qBlocks;
    end
end

nJobs = numel(jobs);
if nJobs == 0
    error('No jobs to run. Fill in DS().blocks for at least one dataset.');
end

% Estimate time
estMin = struct('latulm', 7, 'quasar', 10);  % minutes per block
totalEst = 0;
fprintf('\n  %-5s %-35s %5s %8s %8s\n', 'Job', 'Name', 'Blks', 'Pipeline', 'Est');
fprintf('  %s\n', repmat('-', 1, 70));
for j = 1:nJobs
    nB = numel(jobs(j).blocks);
    e = nB * estMin.(jobs(j).pipeline);
    totalEst = totalEst + e;
    fprintf('  %-5d %-35s %5d %8s %6.0f m\n', j, jobs(j).name, nB, jobs(j).pipeline, e);
end
fprintf('  %s\n', repmat('-', 1, 70));
fprintf('  Pipeline jobs: %d | Post-processing runs after each job\n', nJobs);
fprintf('  Estimated pipeline time: %.0f min (%.1f hours)\n\n', totalEst, totalEst/60);

%% ========================================================================
%  SELECT JOBS
%  ========================================================================
fprintf('  Select: all / 1,3 / 1-%d / q\n', nJobs);
ri = input('  > ', 's'); ri = strtrim(ri);
if strcmpi(ri, 'q'), fprintf('  Cancelled.\n'); return; end
if strcmpi(ri, 'all'), sel = 1:nJobs; else, sel = parse_sel(ri, nJobs); end
if isempty(sel), fprintf('  No valid selection.\n'); return; end

%% ========================================================================
%  GPU CHECK
%  ========================================================================
if C.useGPU
    if gpuDeviceCount > 0
        g = gpuDevice; reset(g);
        fprintf('\n  GPU: %s (%.1f GB VRAM)\n\n', g.Name, g.TotalMemory/1e9);
    else
        fprintf('\n  WARNING: No GPU found. Using CPU.\n\n');
        C.useGPU = false;
    end
end

%% ========================================================================
%  EXECUTE PIPELINE JOBS
%  ========================================================================
batchTimer = tic;

for iJ = 1:numel(sel)
    j = sel(iJ);
    ds = DS(jobs(j).dsIdx);
    jobBlocks = jobs(j).blocks;  % This job's specific block list
    
    fprintf('\n================================================================\n');
    fprintf('  JOB %d/%d: %s\n', iJ, numel(sel), jobs(j).name);
    fprintf('  Pipeline:   %s\n', jobs(j).pipeline);
    fprintf('  Blocks:     %d\n', numel(jobBlocks));
    fprintf('  Output:     %s\n', jobs(j).outDir);
    fprintf('================================================================\n');
    
    if isempty(jobBlocks)
        fprintf('  SKIPPED: no blocks.\n'); continue;
    end
    
    if ~exist(jobs(j).outDir, 'dir'), mkdir(jobs(j).outDir); end
    jobTimer = tic;
    
    try
        cfg = C;
        cfg.runName = jobs(j).name;
        
        % Apply per-dataset config overrides if defined
        if isfield(ds, 'configOverride') && isstruct(ds.configOverride) && ~isempty(fieldnames(ds.configOverride))
            cfg = merge_struct(cfg, ds.configOverride);
            fprintf('  Config overrides applied.\n');
        end
        
        % --- Load metadata from first block ---
        fprintf('  Loading metadata...\n');
        [VT,P,T,~] = VsiVadaDataRead(ds.dataFolder, jobBlocks{1}, ...
            1:cfg.eventsPerFrame, cfg.modeName);
        
        rp = T.ArrayPitch;
        if rp==0||isnan(rp), pm=cfg.pitchOverride;
        elseif rp<10, pm=rp; else, pm=rp/1000; end
        fs=P.SampleFreq; doff=P.DepthOffset;
        cc=P.SoSMedia; if cc==0, cc=1540; end
        if ~isempty(cfg.sosOverride), cc=cfg.sosOverride; end
        
        fprintf('  Probe: %s | pitch=%.3f mm | Fs=%.0f MHz | SoS=%.0f m/s\n', ...
            T.Name, pm, fs, cc);
        
        % Event structure (auto-detect)
        aO = zeros(cfg.eventsPerFrame,1); pO = aO;
        for ev = 1:cfg.eventsPerFrame
            if isfield(VT(ev).TxDelay,'angle'), aO(ev)=VT(ev).TxDelay.angle; end
            if isfield(VT(ev).Waveform,'Channel')&&isfield(VT(ev).Waveform.Channel(1),'invert')
                pO(ev)=VT(ev).Waveform.Channel(1).invert; end
        end
        uA = unique(aO,'stable');
        AP = struct('angle',{},'posIdx',{},'negIdx',{},'rxElements',{});
        for a = 1:numel(uA)
            AP(a).angle=uA(a); ix=find(aO==uA(a));
            if numel(ix)>=2
                if pO(ix(1))==0, AP(a).posIdx=ix(1); AP(a).negIdx=ix(2);
                else, AP(a).posIdx=ix(2); AP(a).negIdx=ix(1); end
            else
                AP(a).posIdx=ix(1); AP(a).negIdx=[];
            end
            AP(a).rxElements=VT(AP(a).posIdx).Elements;
        end
        cfg.numAngles = numel(uA);
        hasPI = numel(ix) >= 2;  % Last angle's count
        cfg.hasPI = hasPI;
        if ~hasPI
            cfg.numPolarities = 1;
            cfg.eventsPerFrame = cfg.numAngles;
        end
        
        zi = find([AP.angle]==0);
        if isempty(zi), [~,zi]=min(abs([AP.angle])); end
        cfg.zeroAngleIdx=zi; cfg.zeroEventIdx=AP(zi).posIdx;
        
        % TX frequency
        txF = 6;
        if isfield(VT(1).Waveform,'Channel')&&~isempty(VT(1).Waveform.Channel)
            txF = VT(1).Waveform.Channel(1).frequency;
        end
        nRx = numel(AP(zi).rxElements);
        clear VT;
        
        % Frame rate
        try
            [VTmp,~,~,~] = VsiVadaDataRead(ds.dataFolder, jobBlocks{1}, ...
                1:(2*cfg.eventsPerFrame), cfg.modeName);
            framePeriod_ms = VTmp(cfg.eventsPerFrame+1).Timestamp - VTmp(1).Timestamp;
            frameRate = 1000 / framePeriod_ms;
            clear VTmp;
        catch
            frameRate = 1000 / (cfg.eventsPerFrame * 0.150);
        end
        fprintf('  Frame rate: %.0f Hz\n', frameRate);
        
        % Grid + delay tables
        eP = ((1:T.ArrayNumElements)-(T.ArrayNumElements+1)/2)*pm;
        xG = cfg.bf.xRange(1):cfg.bf.dx:cfg.bf.xRange(2);
        zG = cfg.bf.zRange(1):cfg.bf.dz:cfg.bf.zRange(2);
        nX = numel(xG); nZ = numel(zG);
        rxP = eP(AP(zi).rxElements);
        
        dT = {beamform_planewave_gpu([], rxP, AP(zi).angle, xG, zG, fs, cc, doff, [])};
        
        % Background
        bgIQ = [];
        if ~isempty(ds.bgFile)
            fprintf('  Computing background...\n');
            try
                [VB,~,~,BC] = VsiVadaDataRead(ds.bgFolder, ds.bgFile, ...
                    1:cfg.eventsPerFrame, cfg.modeName);
                nBE = numel(BC.PulseSequences(1).Events);
                nBF = min(floor(nBE/cfg.eventsPerFrame), 500);
                clear VB;
                [VB,~,~,~] = VsiVadaDataRead(ds.bgFolder, ds.bgFile, ...
                    1:nBF*cfg.eventsPerFrame, cfg.modeName);
                acc = complex(zeros(nZ,nX,'single'));
                for f = 1:nBF
                    ei = (f-1)*cfg.eventsPerFrame + cfg.zeroEventIdx;
                    acc = acc + single(beamform_planewave_gpu(single(VB(ei).Data), ...
                        rxP, AP(zi).angle, xG, zG, fs, cc, doff, dT{1}));
                end
                clear VB; bgIQ = acc/nBF;
                fprintf('    Background: %d frames\n', nBF);
            catch ME
                fprintf('    WARNING: Background failed: %s\n', ME.message);
            end
        end
        
        % --- DISPATCH PIPELINE ---
        switch jobs(j).pipeline
            case 'latulm'
                run_latulm(ds, jobBlocks, cfg, AP, eP, dT, xG, zG, nX, nZ, ...
                    fs, cc, doff, bgIQ, frameRate, jobs(j).outDir);
            case 'quasar'
                run_quasar(ds, jobBlocks, cfg, AP, eP, dT, xG, zG, nX, nZ, ...
                    fs, cc, doff, bgIQ, txF, pm, nRx, frameRate, jobs(j).outDir);
        end
        
    catch ME
        fprintf('\n  FATAL ERROR: %s\n', ME.message);
        for k = 1:min(5, numel(ME.stack))
            fprintf('    at %s (line %d)\n', ME.stack(k).name, ME.stack(k).line);
        end
        fprintf('  Job completed in %.1f min\n', toc(jobTimer)/60);
        if C.useGPU, try reset(gpuDevice); catch, end; end
        continue;  % Skip post-processing for failed jobs
    end
    
    fprintf('  Pipeline completed in %.1f min\n', toc(jobTimer)/60);
    if C.useGPU, try reset(gpuDevice); catch, end; end
    
    % === POST-PROCESSING (runs immediately after each job) ===
    fprintf('\n  --- Post-processing: %s ---\n', jobs(j).name);
    
    % Merge postproc overrides for this dataset
    pp = C.postproc;
    if isfield(ds, 'configOverride') && isstruct(ds.configOverride) && isfield(ds.configOverride, 'postproc')
        pp = merge_struct(pp, ds.configOverride.postproc);
    end
    
    latulmFile = fullfile(ds.outputBase, 'LAT_ULM_Results', 'LAT_ULM_results.mat');
    quasarFile = fullfile(ds.outputBase, 'QUASAR_Results', 'QUASAR_results.mat');
    
    if strcmp(jobs(j).pipeline, 'latulm')
        if ds.runPresentation
            fprintf('    Generating presentation figures...\n');
            try, gen_presentation(latulmFile, ds, pp);
            catch ME, fprintf('      ERROR: %s\n', ME.message); end
        end
        if ds.runAnimation
            fprintf('    Generating animation...\n');
            try, gen_animation(latulmFile, ds, pp);
            catch ME, fprintf('      ERROR: %s\n', ME.message); end
        end
    end
    
    if strcmp(jobs(j).pipeline, 'quasar')
        if ds.runComparison && exist(latulmFile, 'file')
            fprintf('    Generating pipeline comparison...\n');
            try, gen_comparison(latulmFile, quasarFile, ds, pp);
            catch ME, fprintf('      ERROR: %s\n', ME.message); end
        elseif ds.runComparison
            fprintf('    Skipping comparison: LAT-ULM results not yet available.\n');
        end
    end
    
    fprintf('  Job + post-processing done.\n');
end

totalElapsed = toc(batchTimer);
fprintf('\n================================================================\n');
fprintf('  BATCH COMPLETE: %.1f min (%.1f hours)\n', totalElapsed/60, totalElapsed/3600);
fprintf('================================================================\n');

%% ========================================================================
%  LAT-ULM PIPELINE
%  ========================================================================
function run_latulm(ds, jobBlocks, cfg, AP, eP, dT, xG, zG, nX, nZ, fs, cc, doff, bgIQ, frameRate, outDir)
    aL=[]; aT={}; gOff=0; tF=0;
    for iB = 1:numel(jobBlocks)
        tail = jobBlocks{iB}; if numel(tail)>35, tail=['...' tail(end-32:end)]; end
        fprintf('    Block %d/%d: %s\n', iB, numel(jobBlocks), tail);
        
        br = process_single_block(ds.dataFolder, jobBlocks{iB}, cfg.modeName, ...
            cfg, AP, eP, dT, xG, zG, nX, nZ, fs, cc, doff, bgIQ);
        
        if ~isempty(br.localizations)
            br.localizations(:,4) = br.localizations(:,4) + gOff;
            aL = [aL; br.localizations]; %#ok
        end
        for iT = 1:numel(br.tracks)
            t = br.tracks{iT}; t(:,4) = t(:,4) + gOff; aT{end+1} = t; %#ok
        end
        gOff = gOff + br.numFrames; tF = tF + br.numFrames;
    end
    
    R.localizations=aL; R.tracks=aT; R.totalFrames=tF;
    R.config=cfg; R.xGrid=xG; R.zGrid=zG; R.frameRate=frameRate;
    save(fullfile(outDir, 'LAT_ULM_results.mat'), 'R', '-v7.3');
    fprintf('  Saved: %d locs, %d tracks, %d frames\n', size(aL,1), numel(aT), tF);
end

%% ========================================================================
%  QUASAR PIPELINE
%  ========================================================================
function run_quasar(ds, jobBlocks, cfg, AP, eP, dT, xG, zG, nX, nZ, fs, cc, doff, ...
        bgIQ, txF, pm, nRx, frameRate, outDir)
    
    sf=cfg.sushi.srFactor; sp=cfg.bf.dx/sf;
    sNZ=nZ*sf; sNX=nX*sf;
    lambda_mm = (cc*1e-3)/txF;
    [psf,pp] = build_sushi_psf(txF,cc,pm,nRx,sp,[sNZ,sNX],0);
    
    sD=zeros(sNZ,sNX,'single'); qD=sD; nE=0;
    vA=zeros(sNZ,sNX,'single'); vC=vA;
    eS=cfg.ensemble.size; eSt=eS-cfg.ensemble.overlap;
    fO.maxIter=cfg.sushi.maxIter; fO.nonNeg=cfg.sushi.nonNeg; fO.useGPU=cfg.useGPU;
    qO.maxIterCG=cfg.quasar.maxIterCG; qO.supportThresh=cfg.quasar.supportThresh; qO.useGPU=cfg.useGPU;
    zi=cfg.zeroAngleIdx;
    rxP = eP(AP(zi).rxElements);
    
    for iB = 1:numel(jobBlocks)
        bt = tic;
        tail = jobBlocks{iB}; if numel(tail)>35, tail=['...' tail(end-32:end)]; end
        fprintf('    Block %d/%d: %s\n', iB, numel(jobBlocks), tail);
        
        [~,~,~,BC] = VsiVadaDataRead(ds.dataFolder, jobBlocks{iB}, ...
            1:cfg.eventsPerFrame, cfg.modeName);
        nEv = numel(BC.PulseSequences(1).Events);
        nCF = floor(nEv/cfg.eventsPerFrame);
        nBE = floor((nCF-eS)/eSt)+1;
        
        IQ = zeros(nZ,nX,nCF,'single');
        nCh = ceil(nCF/cfg.chunkSize);
        
        fprintf('      Beamforming (%d chunks)...', nCh);
        for iC = 1:nCh
            f1=(iC-1)*cfg.chunkSize+1; f2=min(iC*cfg.chunkSize,nCF); nF=f2-f1+1;
            e1=(f1-1)*cfg.eventsPerFrame+1; e2=f2*cfg.eventsPerFrame;
            try [VC,~,~,~]=VsiVadaDataRead(ds.dataFolder,jobBlocks{iB},e1:e2,cfg.modeName);
            catch, continue; end
            for f=1:nF
                ei=(f-1)*cfg.eventsPerFrame+cfg.zeroEventIdx;
                IQ(:,:,f1+f-1)=single(beamform_planewave_gpu(single(VC(ei).Data),rxP,...
                    AP(zi).angle,xG,zG,fs,cc,doff,dT{1}));
            end
            clear VC;
        end
        fprintf(' done\n');
        
        if ~isempty(bgIQ)
            for f=1:nCF, IQ(:,:,f)=IQ(:,:,f)-bgIQ; end
        end
        
        for iE = 1:nBE
            s1=(iE-1)*eSt+1; s2=s1+eS-1; if s2>nCF, break; end
            ensIQ = svd_clutter_filter_gpu(IQ(:,:,s1:s2), cfg.svd.cutoffLow, ...
                cfg.svd.cutoffHigh, cfg.useGPU);
            
            if cfg.velocity.enable
                [eV,~] = estimate_decorrelation_velocity(ensIQ, ...
                    1000/(cfg.eventsPerFrame*0.150), lambda_mm, sf, cfg.useGPU);
                vm = eV>0 & ~isnan(eV);
                vA(vm) = vA(vm)+eV(vm); vC(vm) = vC(vm)+1;
            end
            
            pw = mean(abs(ensIQ).^2, 3); clear ensIQ;
            pSR = max(imresize(pw,[sNZ,sNX],'bicubic'), 0);
            nF = max(pSR(:))+eps;
            [xf,~] = sushi_sparse_recovery(pSR/nF, psf, cfg.sushi.lambda, cfg.sushi.method, fO);
            sD = sD + xf*nF;
            if cfg.quasar.enable
                [xq,~] = quasar_refit(pSR/nF, psf, xf, qO);
                qD = qD + xq*nF;
            end
            nE = nE+1;
        end
        clear IQ;
        fprintf('      %d ensembles (%.0f sec)\n', min(iE,nBE), toc(bt));
    end
    
    velMap = vA./(vC+eps); velMap(vC==0) = 0;
    srX = linspace(cfg.bf.xRange(1),cfg.bf.xRange(2),sNX);
    srZ = linspace(cfg.bf.zRange(1),cfg.bf.zRange(2),sNZ);
    
    R.sushiDensity=sD; R.quasarDensity=qD; R.velocityMap=velMap;
    R.nEnsembles=nE; R.config=cfg; R.psfParams=pp;
    R.frameRate=1000/(cfg.eventsPerFrame*0.150); R.srX=srX; R.srZ=srZ;
    save(fullfile(outDir,'QUASAR_results.mat'),'R','-v7.3');
    fprintf('  Saved: %d ensembles\n', nE);
end

%% ========================================================================
%  POST-PROCESSING: PRESENTATION FIGURES
%  ========================================================================
function gen_presentation(resultsFile, ds, pp)
    tmp = load(resultsFile);
    if isfield(tmp,'R'), R=tmp.R; elseif isfield(tmp,'results'), R=tmp.results;
    else, fn=fieldnames(tmp); R=tmp.(fn{1}); end
    
    tracks = R.tracks; frameRate = R.frameRate;
    outDir = fileparts(resultsFile);
    
    % Filter slow tracks
    keepMask = true(numel(tracks), 1);
    for i = 1:numel(tracks)
        t = tracks{i};
        if size(t,1) > 1
            dx=diff(t(:,1)); dz=diff(t(:,2)); df=diff(t(:,4));
            dt=df/frameRate; sp=sqrt(dx.^2+dz.^2)./(dt+eps);
            if median(sp) < pp.minSpeed_mm_s, keepMask(i) = false; end
        else
            keepMask(i) = false;
        end
    end
    tracks = tracks(keepMask);
    nTracks = numel(tracks);
    fprintf('    %d tracks after speed filter (>%.1f mm/s)\n', nTracks, pp.minSpeed_mm_s);
    
    if nTracks == 0
        fprintf('    No tracks after filtering — skipping figures.\n');
        return;
    end
    
    allPts = cell2mat(tracks(:));
    xLim = [min(allPts(:,1))-0.5, max(allPts(:,1))+0.5];
    zLim = [min(allPts(:,2))-0.5, max(allPts(:,2))+0.5];
    px = pp.srPixel_um / 1000;
    xEdges = xLim(1):px:xLim(2); zEdges = zLim(1):px:zLim(2);
    nXsr = numel(xEdges)-1; nZsr = numel(zEdges)-1;
    xC = xEdges(1:end-1)+px/2; zC = zEdges(1:end-1)+px/2;
    
    % Build density + velocity maps
    density = zeros(nZsr, nXsr);
    velAccum = zeros(nZsr, nXsr); velWeight = zeros(nZsr, nXsr);
    dirAccum = zeros(nZsr, nXsr);
    medSpeeds = zeros(nTracks, 1); netDirs = zeros(nTracks, 1);
    
    for i = 1:nTracks
        t = tracks{i};
        dx=diff(t(:,1)); dz=diff(t(:,2)); df=diff(t(:,4));
        dt=df/frameRate; sp=sqrt(dx.^2+dz.^2)./(dt+eps);
        medSpeeds(i) = median(sp);
        netDirs(i) = sign(t(end,2)-t(1,2));
        
        for k = 1:size(t,1)
            xi = find(xEdges(1:end-1)<=t(k,1), 1, 'last');
            zi = find(zEdges(1:end-1)<=t(k,2), 1, 'last');
            if ~isempty(xi) && ~isempty(zi) && xi<=nXsr && zi<=nZsr
                density(zi,xi) = density(zi,xi) + 1;
                velAccum(zi,xi) = velAccum(zi,xi) + medSpeeds(i);
                velWeight(zi,xi) = velWeight(zi,xi) + 1;
                dirAccum(zi,xi) = dirAccum(zi,xi) + netDirs(i);
            end
        end
    end
    
    velMap = velAccum ./ max(velWeight,1); velMap(velWeight==0)=NaN;
    kern = fspecial('gaussian',[5 5],1.0);
    velSmooth = nanconv_local(velMap, kern);
    velSmooth(velWeight==0) = NaN;
    velDir = velSmooth .* sign(dirAccum); velDir(velWeight==0) = NaN;
    
    % Three-panel
    fig = figure('Visible','off','Position',[50 100 1800 550],'Color','w');
    
    ax1=subplot(1,3,1);
    imagesc(xC,zC,log10(density+1)); axis image; colormap(ax1,hot); cb1=colorbar;
    ylabel(cb1,'log_{10}(count+1)','Color','k');
    xlabel('Lateral [mm]','Color','k'); ylabel('Axial [mm]','Color','k');
    title(sprintf('%s\nDensity (%d tracks)', ds.name, nTracks),'Color','k');
    set(ax1,'FontSize',10,'XColor','k','YColor','k');
    
    ax2=subplot(1,3,2); set(ax2,'Color','k');
    halfN=128;
    b2w=[linspace(0.2,1,halfN)',linspace(0.3,1,halfN)',ones(halfN,1)];
    w2r=[ones(halfN,1),linspace(1,0.2,halfN)',linspace(1,0.2,halfN)'];
    hImg=imagesc(xC,zC,velDir); set(hImg,'AlphaData',~isnan(velDir));
    axis image; colormap(ax2,[b2w;w2r]); cb2=colorbar; caxis([-pp.velMax pp.velMax]);
    ylabel(cb2,'Speed [mm/s] (blue=up, red=down)','Color','k');
    xlabel('Lateral [mm]','Color','k'); ylabel('Axial [mm]','Color','k');
    velAnnot='';
    for i=1:numel(ds.tubes)
        area=pi*(ds.tubes(i).ID_mm/2)^2;
        Q=ds.tubes(i).rate_mL_min*1000/60;
        vExp=Q/area;
        velAnnot=[velAnnot,sprintf('%s: %.1f mm/s  ',ds.tubes(i).name,vExp)]; %#ok
    end
    title(sprintf('Velocity\n%s',velAnnot),'Color','k');
    set(ax2,'FontSize',10,'XColor','k','YColor','k');
    
    ax3=subplot(1,3,3); hold on;
    cm=lines(min(nTracks,256));
    for i=1:nTracks
        t=tracks{i}; col=cm(mod(i-1,size(cm,1))+1,:);
        plot(t(:,1),t(:,2),'-','Color',[col 0.7],'LineWidth',0.8);
    end
    set(gca,'YDir','reverse'); axis equal; xlim(xLim); ylim(zLim);
    xlabel('Lateral [mm]','Color','k'); ylabel('Axial [mm]','Color','k');
    title(sprintf('Tracks (%d)',nTracks),'Color','k'); grid on;
    set(ax3,'FontSize',10,'XColor','k','YColor','k','GridAlpha',0.15);
    
    exportgraphics(fig,fullfile(outDir,'three_panel.png'),'Resolution',300);
    close(fig);
    fprintf('    Saved: three_panel.png\n');
end

%% ========================================================================
%  POST-PROCESSING: ANIMATION
%  ========================================================================
function gen_animation(resultsFile, ds, pp)
    tmp = load(resultsFile);
    if isfield(tmp,'R'), R=tmp.R; elseif isfield(tmp,'results'), R=tmp.results;
    else, fn=fieldnames(tmp); R=tmp.(fn{1}); end
    
    tracks = R.tracks; frameRate = R.frameRate;
    outDir = fileparts(resultsFile);
    nTracks = numel(tracks);
    
    % Speed filter
    keepMask = true(nTracks,1);
    medSpeeds = zeros(nTracks,1); netDirs = zeros(nTracks,1);
    for i=1:nTracks
        t=tracks{i};
        if size(t,1)>1
            dx=diff(t(:,1));dz=diff(t(:,2));df=diff(t(:,4));
            dt=df/frameRate; sp=sqrt(dx.^2+dz.^2)./(dt+eps);
            medSpeeds(i)=median(sp); netDirs(i)=sign(t(end,2)-t(1,2));
            if medSpeeds(i)<pp.minSpeed_mm_s, keepMask(i)=false; end
        else, keepMask(i)=false; end
    end
    tracks=tracks(keepMask); medSpeeds=medSpeeds(keepMask); netDirs=netDirs(keepMask);
    nTracks=numel(tracks);
    
    if nTracks == 0
        fprintf('    No tracks after filtering — skipping animation.\n');
        return;
    end
    
    % Smooth
    w = pp.anim_smoothWindow;
    for i=1:nTracks
        t=tracks{i};
        if size(t,1)>w
            tracks{i}(:,1)=movmean(t(:,1),w);
            tracks{i}(:,2)=movmean(t(:,2),w);
        end
    end
    
    % Frame ranges
    tRanges = zeros(nTracks,2);
    for i=1:nTracks
        tRanges(i,:) = [min(tracks{i}(:,4)), max(tracks{i}(:,4))];
    end
    allPts = cell2mat(tracks(:));
    xLim = [min(allPts(:,1))-0.5, max(allPts(:,1))+0.5];
    zLim = [min(allPts(:,2))-0.5, max(allPts(:,2))+0.5];
    
    fMin=min(tRanges(:,1)); fMax=max(tRanges(:,2));
    step=3; maxDF=pp.anim_duration_s*pp.anim_fps*step;
    if (fMax-fMin)>maxDF, fMax=fMin+maxDF; end
    aFrames = fMin:step:fMax;
    speedMax = prctile(medSpeeds,95); if speedMax==0, speedMax=pp.velMax; end
    
    fig=figure('Position',[100 100 800 900],'Color','k','MenuBar','none','Visible','off');
    ax=axes('Parent',fig,'Color','k','XColor','w','YColor','w');
    hold(ax,'on'); set(ax,'YDir','reverse');
    xlim(ax,xLim); ylim(ax,zLim); axis(ax,'equal');
    xlabel(ax,'Lateral [mm]','Color','w'); ylabel(ax,'Axial [mm]','Color','w');
    
    vidFile = fullfile(outDir, sprintf('%s_animation.mp4', ds.shortName));
    v=VideoWriter(vidFile,'MPEG-4'); v.FrameRate=pp.anim_fps; v.Quality=95; open(v);
    
    for iA=1:numel(aFrames)
        cf=aFrames(iA); cla(ax);
        ws=cf-40;
        active = find((tRanges(:,1)<=cf) & (tRanges(:,2)>=ws));
        for ii=1:numel(active)
            ti=active(ii); t=tracks{ti};
            inW=(t(:,4)>=ws)&(t(:,4)<=cf); tw=t(inW,:);
            if isempty(tw), continue; end
            frac=min(medSpeeds(ti)/speedMax,1);
            col=[frac,0.1,1-frac];
            if size(tw,1)>1, plot(ax,tw(:,1),tw(:,2),'-','Color',[col 0.5],'LineWidth',1.5); end
            [~,li]=max(tw(:,4)); plot(ax,tw(li,1),tw(li,2),'.','Color',col,'MarkerSize',10);
        end
        title(ax,sprintf('%s | t=%.2fs',ds.name,(cf-fMin)/frameRate),'Color','w');
        drawnow limitrate; writeVideo(v,getframe(fig));
    end
    close(v); close(fig);
    fprintf('    Saved: %s (%.1f s)\n', ds.shortName, numel(aFrames)/pp.anim_fps);
end

%% ========================================================================
%  POST-PROCESSING: LAT-ULM vs QUASAR COMPARISON
%  ========================================================================
function gen_comparison(latulmFile, quasarFile, ds, pp)
    % Load LAT-ULM
    tmp = load(latulmFile);
    if isfield(tmp,'R'), LR=tmp.R; else, fn=fieldnames(tmp); LR=tmp.(fn{1}); end
    
    % Load QUASAR
    tmp = load(quasarFile);
    if isfield(tmp,'R'), QR=tmp.R; else, fn=fieldnames(tmp); QR=tmp.(fn{1}); end
    
    outDir = ds.outputBase;
    if ~exist(outDir,'dir'), mkdir(outDir); end
    
    tracks = LR.tracks; frameRate = LR.frameRate;
    
    % Filter LAT-ULM tracks
    keepMask = true(numel(tracks),1);
    for i=1:numel(tracks)
        t=tracks{i};
        if size(t,1)>1
            dx=diff(t(:,1));dz=diff(t(:,2));df=diff(t(:,4));
            dt=df/frameRate; sp=sqrt(dx.^2+dz.^2)./(dt+eps);
            if median(sp)<pp.minSpeed_mm_s, keepMask(i)=false; end
        else, keepMask(i)=false; end
    end
    tracks=tracks(keepMask);
    
    if isempty(tracks)
        fprintf('    No LAT-ULM tracks after filtering — skipping comparison.\n');
        return;
    end
    
    % Build LAT-ULM density on same grid as QUASAR
    srX = QR.srX; srZ = QR.srZ;
    px = srX(2)-srX(1); nXsr=numel(srX); nZsr=numel(srZ);
    xEdges = [srX - px/2, srX(end)+px/2];
    zEdges = [srZ - px/2, srZ(end)+px/2];
    
    latulmDens = zeros(nZsr, nXsr);
    for i=1:numel(tracks)
        t=tracks{i};
        for k=1:size(t,1)
            xi=find(xEdges(1:end-1)<=t(k,1),1,'last');
            zi=find(zEdges(1:end-1)<=t(k,2),1,'last');
            if ~isempty(xi)&&~isempty(zi)&&xi<=nXsr&&zi<=nZsr
                latulmDens(zi,xi)=latulmDens(zi,xi)+1;
            end
        end
    end
    
    % --- Figure: 2x2 comparison ---
    fig = figure('Visible','off','Position',[50 50 1200 900],'Color','w');
    
    % LAT-ULM density
    ax=subplot(2,2,1);
    imagesc(srX, srZ, log10(latulmDens+1));
    axis image; colormap(gca, hot); cb=colorbar;
    ylabel(cb,'log_{10}(count+1)','Color','k');
    title(sprintf('LAT-ULM Density\n%d tracks', numel(tracks)),'Color','k');
    xlabel('Lateral [mm]','Color','k'); ylabel('Axial [mm]','Color','k');
    set(ax,'FontSize',10,'XColor','k','YColor','k');
    
    % QUASAR density
    ax=subplot(2,2,2);
    imagesc(srX, srZ, log10(QR.quasarDensity+1));
    axis image; colormap(gca, hot); cb=colorbar;
    ylabel(cb,'log_{10}(count+1)','Color','k');
    title(sprintf('QUASAR Density\n%d ensembles', QR.nEnsembles),'Color','k');
    xlabel('Lateral [mm]','Color','k'); ylabel('Axial [mm]','Color','k');
    set(ax,'FontSize',10,'XColor','k','YColor','k');
    
    % LAT-ULM velocity (track-painted)
    ax=subplot(2,2,3);
    velAccum=zeros(nZsr,nXsr); velW=zeros(nZsr,nXsr);
    for i=1:numel(tracks)
        t=tracks{i}; if size(t,1)<2, continue; end
        dx=diff(t(:,1));dz=diff(t(:,2));df=diff(t(:,4));
        dt=df/frameRate; sp=sqrt(dx.^2+dz.^2)./(dt+eps);
        ms=median(sp);
        for k=1:size(t,1)
            xi=find(xEdges(1:end-1)<=t(k,1),1,'last');
            zi=find(zEdges(1:end-1)<=t(k,2),1,'last');
            if ~isempty(xi)&&~isempty(zi)&&xi<=nXsr&&zi<=nZsr
                velAccum(zi,xi)=velAccum(zi,xi)+ms; velW(zi,xi)=velW(zi,xi)+1;
            end
        end
    end
    vMap=velAccum./max(velW,1); vMap(velW==0)=NaN;
    set(ax,'Color','k');
    hImg=imagesc(srX,srZ,vMap); set(hImg,'AlphaData',~isnan(vMap));
    axis image; cb=colorbar; caxis([0 pp.velMax]); colormap(gca,jet);
    ylabel(cb,'mm/s','Color','k');
    title('LAT-ULM Velocity [mm/s]','Color','k');
    xlabel('Lateral [mm]','Color','k'); ylabel('Axial [mm]','Color','k');
    set(ax,'FontSize',10,'XColor','k','YColor','k');
    
    % QUASAR velocity
    ax=subplot(2,2,4);
    qVel = QR.velocityMap; qVel(qVel==0) = NaN;
    set(ax,'Color','k');
    hImg=imagesc(srX, srZ, qVel); set(hImg,'AlphaData',~isnan(qVel));
    axis image; cb=colorbar; caxis([0 pp.velMax]); colormap(gca,jet);
    ylabel(cb,'mm/s','Color','k');
    title('QUASAR Velocity [mm/s]','Color','k');
    xlabel('Lateral [mm]','Color','k'); ylabel('Axial [mm]','Color','k');
    set(ax,'FontSize',10,'XColor','k','YColor','k');
    
    sgtitle(sprintf('%s: LAT-ULM vs QUASAR', ds.name),'Color','k');
    exportgraphics(fig,fullfile(outDir,'pipeline_comparison.png'),'Resolution',300);
    close(fig);
    
    % --- Lateral profile comparison ---
    fig2 = figure('Visible','off','Position',[100 100 800 400],'Color','w');
    depth = pp.velProfile_depth_mm;
    
    zMaskL = abs(srZ-depth) <= 0.25;
    profL = mean(latulmDens(zMaskL,:),1);
    profL = profL / (max(profL)+eps);
    profQ = mean(QR.quasarDensity(zMaskL,:),1);
    profQ = profQ / (max(profQ)+eps);
    
    plot(srX, profL, 'b-', 'LineWidth', 1.5, 'DisplayName', 'LAT-ULM');
    hold on;
    plot(srX, profQ, 'r-', 'LineWidth', 1.5, 'DisplayName', 'QUASAR');
    yline(0.5, 'k--', 'FWHM');
    xlabel('Lateral [mm]','Color','k','FontSize',12);
    ylabel('Normalized density','Color','k','FontSize',12);
    title(sprintf('Lateral Profile at z = %.1f mm', depth),'Color','k','FontSize',13);
    legend('Location','best','TextColor','k'); grid on;
    set(gca,'FontSize',11,'XColor','k','YColor','k','Box','on','GridAlpha',0.15);
    
    exportgraphics(fig2,fullfile(outDir,'profile_comparison.png'),'Resolution',300);
    close(fig2);
    
    fprintf('    Saved: pipeline_comparison.png, profile_comparison.png\n');
end

%% ========================================================================
%  UTILITIES
%  ========================================================================
function out = merge_struct(base, override)
% MERGE_STRUCT  Recursively merge override fields into base struct.
    out = base;
    fns = fieldnames(override);
    for i = 1:numel(fns)
        f = fns{i};
        if isstruct(override.(f)) && isfield(base, f) && isstruct(base.(f))
            out.(f) = merge_struct(base.(f), override.(f));
        else
            out.(f) = override.(f);
        end
    end
end

function sel = parse_sel(str, nMax)
    sel = [];
    parts = strsplit(strtrim(str), ',');
    for i = 1:numel(parts)
        p = strtrim(parts{i});
        if contains(p, '-')
            rng = sscanf(p, '%d-%d');
            if numel(rng)==2, sel = [sel, rng(1):rng(2)]; end %#ok
        else
            v = str2double(p);
            if ~isnan(v), sel = [sel, v]; end %#ok
        end
    end
    sel = unique(sel(sel>=1 & sel<=nMax));
end

function result = nanconv_local(A, kernel)
    A_zero = A; A_zero(isnan(A)) = 0;
    weight = double(~isnan(A));
    result = conv2(A_zero, kernel, 'same') ./ max(conv2(weight, kernel, 'same'), eps);
    result(isnan(A)) = NaN;
end