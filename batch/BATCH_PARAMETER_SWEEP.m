%% BATCH_PARAMETER_SWEEP.m
% Multi-dataset parameter sweep: beamforms ALL blocks ONCE per dataset,
% caches pre-SVD IQ to disk, then sweeps tracking parameters with
% three-panel output for every combination. Designed for overnight runs.
%
% DATASETS: Primary (higher conc), Primary_Low (5e5), Mac+MB
% PIPELINES: LAT-ULM (tracking sweep) + QUASAR (per SVD cutoff)
%
% OUTPUT STRUCTURE (per dataset):
%   <outputBase>/
%     IQ_cache/                     ← cached pre-SVD IQ (one .mat per block)
%     svd2_thr5_len10_gap5_.../     ← one folder per param combo
%       LAT_ULM_results.mat
%       three_panel.png
%     QUASAR_svd2/                  ← QUASAR results per SVD cutoff
%       QUASAR_results.mat
%     sweep_summary.mat             ← all metrics for comparison
%     sweep_comparison.png          ← marginalized heatmaps
%
% Author: Eli Wirth-Apley / Sun Lab, Northeastern University, March 2026

clearvars; close all; clc;

%% ========================================================================
%  CONFIGURATION
%  ========================================================================
baseDir         = 'C:\path\to\ULM3';
macmbDir        = 'C:\path\to\data\Series1';
vadaScriptsPath = 'C:\path\to\VevoF2\RFExampleScriptsF2\Matlab';
addpath(genpath(vadaScriptsPath));
pipelineDir = fileparts(mfilename('fullpath'));
if isempty(pipelineDir), pipelineDir = pwd; end
addpath(pipelineDir);

% --- Datasets ---
DS = struct(); nDS = 0;

% Dataset 1: Primary (higher concentration, ~2e6)
nDS = nDS+1;
DS(nDS).name       = 'Primary (Higher Conc)';
DS(nDS).shortName  = 'primary';
DS(nDS).dataFolder = fullfile(baseDir, 'Primary');
DS(nDS).bgFile     = '';  % disabled — SVD handles clutter, bg sub can inject artifacts
DS(nDS).bgFolder   = '';
DS(nDS).blocks     = { ...
    '2026-02-25-23-28-57_2026-02-25-22-43-32', ...
    '2026-02-25-23-28-57_2026-02-25-22-45-05', ...
    '2026-02-25-23-28-57_2026-02-25-22-46-34', ...
    '2026-02-25-23-28-57_2026-02-25-22-48-13', ...
    '2026-02-25-23-28-57_2026-02-25-22-49-42', ...
    '2026-02-25-23-28-57_2026-02-25-22-51-16', ...
    '2026-02-25-23-28-57_2026-02-25-22-52-48', ...
    '2026-02-25-23-28-57_2026-02-25-22-54-28', ...
    '2026-02-25-23-28-57_2026-02-25-22-55-58'};
DS(nDS).quasarBlocks = { ...
    '2026-02-25-23-28-57_2026-02-25-22-45-05', ...
    '2026-02-25-23-28-57_2026-02-25-22-46-34'};
DS(nDS).outputBase = fullfile(baseDir, 'Primary', 'Sweep_Results');
DS(nDS).runLATULM  = true;
DS(nDS).runQUASAR  = true;
DS(nDS).tubes = struct( ...
    'name',       {'Small tube', 'Large tube'}, ...
    'ID_mm',      {0.31,          0.51}, ...
    'rate_mL_min',{0.08,          0.16}, ...
    'approxLateral_mm', {-0.6,    0.5});
DS(nDS).overrides = struct();  % No overrides — use defaults

% Dataset 2: Primary_Low (5e5 concentration)
nDS = nDS+1;
DS(nDS).name       = 'Phantom 5e5 Initial Speed';
DS(nDS).shortName  = 'primary_low';
DS(nDS).dataFolder = fullfile(baseDir, 'Primary_Low');
DS(nDS).bgFile     = '';  % disabled
DS(nDS).bgFolder   = '';
DS(nDS).blocks     = { ...
    '2026-03-10-10-56-09_2026-02-25-22-16-10', ...
    '2026-03-10-10-56-09_2026-02-25-22-18-22', ...
    '2026-03-10-10-56-09_2026-02-25-22-19-50', ...
    '2026-03-10-10-56-09_2026-02-25-22-21-30', ...
    '2026-03-10-10-56-09_2026-02-25-22-22-56', ...
    '2026-03-10-10-56-09_2026-02-25-22-24-53', ...
    '2026-03-10-10-56-09_2026-02-25-22-26-24', ...
    '2026-03-10-10-56-09_2026-02-25-22-28-06', ...
    '2026-03-10-10-56-09_2026-02-25-22-29-54', ...
    '2026-03-10-10-56-09_2026-02-25-22-31-28'};
DS(nDS).quasarBlocks = { ...
    '2026-03-10-10-56-09_2026-02-25-22-19-50', ...
    '2026-03-10-10-56-09_2026-02-25-22-21-30'};
DS(nDS).outputBase = fullfile(baseDir, 'Primary_Low', 'Sweep_Results');
DS(nDS).runLATULM  = true;
DS(nDS).runQUASAR  = true;
DS(nDS).tubes = struct( ...
    'name',       {'Small tube', 'Large tube'}, ...
    'ID_mm',      {0.31,          0.51}, ...
    'rate_mL_min',{0.08,          0.16}, ...
    'approxLateral_mm', {-0.6,    0.5});
DS(nDS).overrides = struct();  % No overrides

% Dataset 3: Mac+MB (macrophage-loaded MBs, all blocks)
nDS = nDS+1;
DS(nDS).name       = 'Mac+MB';
DS(nDS).shortName  = 'mac_mb';
DS(nDS).dataFolder = macmbDir;
DS(nDS).bgFile     = '';  % disabled
DS(nDS).bgFolder   = '';
DS(nDS).blocks     = { ...
    'block_1', ...
    'block_2', ...
    'block_3', ...
    'block_4', ...
    'block_5', ...
    'block_6', ...
    'block_7', ...
    'block_8', ...
    'block_9', ...
    'block_10'};
DS(nDS).quasarBlocks = { ...
    'block_2', ...
    'block_3'};
DS(nDS).outputBase = fullfile(macmbDir, 'Sweep_Results');
DS(nDS).runLATULM  = true;
DS(nDS).runQUASAR  = true;
DS(nDS).tubes = struct( ...
    'name',       {'Small tube', 'Large tube'}, ...
    'ID_mm',      {0.31,          0.51}, ...
    'rate_mL_min',{0.08,          0.16}, ...
    'approxLateral_mm', {-0.25,   0.15});
% Mac+MB overrides: coarser resolution, higher threshold, larger gap
DS(nDS).overrides.det.roiSize_px       = 11;
DS(nDS).overrides.det.minSep_mm        = 0.200;
DS(nDS).overrides.sr.pixelSize_um      = 25;
DS(nDS).overrides.postproc.srPixel_um  = 50;
DS(nDS).overrides.postproc.minSpeed_mm_s = 2.0;

% --- Parameter sweep ranges ---
sweep.svdCutoffs      = [2, 3, 5, 8, 12];
sweep.detThresholds   = [3, 5, 8, 12];           % re-added thr=3 (should work with zeroOnly=true)
sweep.minTrackLengths = [5, 10, 15, 20];
sweep.maxGapFrames    = [3, 5, 8];
sweep.maxDisps        = [0.05, 0.07, 0.10];
sweep.processNoises   = [0.002, 0.01];           % trimmed: 0.0005 showed minimal effect
sweep.measNoises      = [0.080, 0.150];           % trimmed: 0.040 showed minimal effect

% --- Fixed parameters (defaults, overridden per-dataset) ---
defaults.modeName       = '.vada';
defaults.useGPU         = true;
defaults.zeroOnly       = true;   % true: avoids grating lobes at pitch/λ≈5.7
defaults.pitchOverride  = [];
defaults.sosOverride    = [];
defaults.chunkSize      = 400;
defaults.eventsPerFrame = [];   % [] = auto-detect
defaults.bf.xRange      = [-2, 2];   % tighter FOV — tubes within ±1mm, margin for grating lobe check
defaults.bf.zRange      = [14, 22];   % tubes at z≈15-18mm
defaults.bf.dx          = 0.025;
defaults.bf.dz          = 0.010;
defaults.det.method     = 'threshold';
defaults.det.minSep_mm  = 0.100;
defaults.det.roiSize_px = 7;
defaults.track.maxDisp_mm           = 0.070;
defaults.track.maxGapFrames         = 5;
defaults.track.kalman.processNoise  = 0.0005;
defaults.track.kalman.measNoise     = 0.080;
defaults.sr.pixelSize_um = 5;
defaults.postproc.minSpeed_mm_s     = 1.5;
defaults.postproc.srPixel_um        = 25;
defaults.postproc.velMax             = 30;
defaults.postproc.fwhm_depths       = [16.0, 16.5, 17.0, 17.5, 18.0];
defaults.postproc.velProfile_depth_mm = 17.0;

% QUASAR parameters
defaults.quasar.lambda   = [0.1, 0.5, 1.0];  % Sweep lambda values
defaults.quasar.srPixel_mm = 0.005;
defaults.quasar.ensembleSize = 200;

%% ========================================================================
%  GPU CHECK
%  ========================================================================
useGPU = defaults.useGPU;
if useGPU
    if gpuDeviceCount > 0
        g = gpuDevice; reset(g);
        fprintf('GPU: %s (%.1f GB VRAM)\n\n', g.Name, g.TotalMemory/1e9);
    else
        fprintf('WARNING: No GPU found. Using CPU.\n\n');
        useGPU = false;
    end
end
defaults.useGPU = useGPU;

nCombos = numel(sweep.svdCutoffs) * numel(sweep.detThresholds) * ...
    numel(sweep.minTrackLengths) * numel(sweep.maxGapFrames) * ...
    numel(sweep.maxDisps) * numel(sweep.processNoises) * numel(sweep.measNoises);

fprintf('================================================================\n');
fprintf('  BATCH PARAMETER SWEEP — %d datasets, %d combos each\n', nDS, nCombos);
fprintf('  SVD: %s | Thr: %s | Len: %s\n', mat2str(sweep.svdCutoffs), ...
    mat2str(sweep.detThresholds), mat2str(sweep.minTrackLengths));
fprintf('  Gap: %s | Disp: %s | ProcNoise: %s | MeasNoise: %s\n', ...
    mat2str(sweep.maxGapFrames), mat2str(sweep.maxDisps), ...
    mat2str(sweep.processNoises), mat2str(sweep.measNoises));
for d = 1:nDS
    fprintf('  [%d] %s — %d blocks, QUASAR: %s\n', d, DS(d).name, ...
        numel(DS(d).blocks), mat2str(DS(d).runQUASAR));
end
fprintf('================================================================\n\n');

masterTimer = tic;

%% ========================================================================
%  MAIN DATASET LOOP
%  ========================================================================
for iDS = 1:nDS

dsTimer = tic;
ds = DS(iDS);
fprintf('\n################################################################\n');
fprintf('  DATASET %d/%d: %s\n', iDS, nDS, ds.name);
fprintf('################################################################\n');

% Merge defaults with per-dataset overrides
fixed = defaults;
if ~isempty(fieldnames(ds.overrides))
    fixed = merge_struct(fixed, ds.overrides);
end

dataFolder = ds.dataFolder;
blocks     = ds.blocks;
bgFile     = ds.bgFile;
bgFolder   = ds.bgFolder;
outputBase = ds.outputBase;
label      = ds.name;
tubes      = ds.tubes;

if ~exist(outputBase, 'dir'), mkdir(outputBase); end
cacheDir = fullfile(outputBase, 'IQ_cache');
if ~exist(cacheDir, 'dir'), mkdir(cacheDir); end

%% ========================================================================
%  LOAD METADATA + DELAY TABLES
%  ========================================================================
fprintf('[1/5] Loading metadata...\n');

% Auto-detect events per frame if not specified
if isempty(fixed.eventsPerFrame)
    numProbe = 30;
    [VadaProbe,~,~,~] = VsiVadaDataRead(dataFolder, blocks{1}, ...
        1:numProbe, fixed.modeName);
    probeAngles = zeros(numProbe,1); probePolar = zeros(numProbe,1);
    for ev = 1:numProbe
        if isfield(VadaProbe(ev).TxDelay,'angle')
            probeAngles(ev) = VadaProbe(ev).TxDelay.angle; end
        if isfield(VadaProbe(ev).Waveform,'Channel') && ...
                isfield(VadaProbe(ev).Waveform.Channel(1),'invert')
            probePolar(ev) = VadaProbe(ev).Waveform.Channel(1).invert; end
    end
    probeSig = probeAngles*10 + probePolar;
    for epf = 1:15
        if mod(numProbe, epf) == 0 || numProbe >= 2*epf
            pat = probeSig(1:epf);
            nRep = floor(numProbe/epf);
            match = true;
            for r = 2:nRep
                if ~isequal(probeSig((r-1)*epf+1:r*epf), pat)
                    match = false; break;
                end
            end
            if match && nRep >= 2
                fixed.eventsPerFrame = epf; break;
            end
        end
    end
    if isempty(fixed.eventsPerFrame), fixed.eventsPerFrame = 6; end
    clear VadaProbe;
    fprintf('  Auto-detected eventsPerFrame = %d\n', fixed.eventsPerFrame);
end

[VT,P,T,~] = VsiVadaDataRead(dataFolder, blocks{1}, ...
    1:fixed.eventsPerFrame, fixed.modeName);

rp = T.ArrayPitch;
if rp==0||isnan(rp)
    if ~isempty(fixed.pitchOverride), pm=fixed.pitchOverride;
    else, pm=0.300; fprintf('  WARNING: ArrayPitch=0, using 0.300 mm default\n');
    end
elseif rp<10, pm=rp; else, pm=rp/1000; end
fs=P.SampleFreq; doff=P.DepthOffset;
cc=P.SoSMedia; if cc==0, cc=1540; end
if ~isempty(fixed.sosOverride), cc=fixed.sosOverride; end

% Event structure
aO = zeros(fixed.eventsPerFrame,1); pO = aO;
for ev = 1:fixed.eventsPerFrame
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
numAngles = numel(AP);
hasPulseInversion = ~isempty(AP(1).negIdx);
zi = find([AP.angle]==0);
if isempty(zi), [~,zi]=min(abs([AP.angle])); end
zeroEventIdx = AP(zi).posIdx;

% Frame rate
try
    [VTmp,~,~,~] = VsiVadaDataRead(dataFolder, blocks{1}, ...
        1:(2*fixed.eventsPerFrame), fixed.modeName);
    framePeriod_ms = VTmp(fixed.eventsPerFrame+1).Timestamp - VTmp(1).Timestamp;
    frameRate = 1000 / framePeriod_ms;
    clear VTmp;
catch
    frameRate = 1000 / (fixed.eventsPerFrame * 0.150);
end
clear VT;

fprintf('  %s: angles=%d, PI=%s, EPF=%d, FR=%.0f Hz\n', ...
    ds.name, numAngles, string(hasPulseInversion), fixed.eventsPerFrame, frameRate);

% Grid
eP = ((1:T.ArrayNumElements)-(T.ArrayNumElements+1)/2)*pm;
xG = fixed.bf.xRange(1):fixed.bf.dx:fixed.bf.xRange(2);
zG = fixed.bf.zRange(1):fixed.bf.dz:fixed.bf.zRange(2);
nX = numel(xG); nZ = numel(zG);
rxP = eP(AP(zi).rxElements);

% Delay table (zero angle only for zeroOnly=true)
if fixed.zeroOnly
    dT = {beamform_planewave_gpu([], rxP, AP(zi).angle, xG, zG, fs, cc, doff, [])};
else
    dT = cell(1, numel(AP));
    for a = 1:numel(AP)
        rxPa = eP(AP(a).rxElements);
        dT{a} = beamform_planewave_gpu([], rxPa, AP(a).angle, xG, zG, fs, cc, doff, []);
    end
end

fprintf('  Probe: %s | pitch=%.3f mm | Fs=%.0f MHz | FR=%.0f Hz\n', T.Name, pm, fs, frameRate);
fprintf('  Grid: %d x %d\n', nX, nZ);

%% ========================================================================
%  BEAMFORM + CACHE (per-block, pre-SVD IQ to disk)
%  ========================================================================
bfTimer = tic;
fprintf('\n[2/5] Beamforming + caching %d blocks...\n', numel(blocks));

% Background
bgIQ = [];
if ~isempty(bgFile)
    bgCacheFile = fullfile(cacheDir, 'bgIQ.mat');
    if exist(bgCacheFile, 'file')
        fprintf('  Loading cached background...\n');
        tmp = load(bgCacheFile, 'bgIQ'); bgIQ = tmp.bgIQ;
    else
        fprintf('  Computing background...\n');
        try
            [VB,~,~,BC] = VsiVadaDataRead(bgFolder, bgFile, ...
                1:fixed.eventsPerFrame, fixed.modeName);
            nBE = numel(BC.PulseSequences(1).Events);
            nBF = min(floor(nBE/fixed.eventsPerFrame), 500);
            clear VB;
            [VB,~,~,~] = VsiVadaDataRead(bgFolder, bgFile, ...
                1:nBF*fixed.eventsPerFrame, fixed.modeName);
            acc = complex(zeros(nZ,nX,'single'));
            for f = 1:nBF
                baseEv = (f-1)*fixed.eventsPerFrame;
                if fixed.zeroOnly
                    ei = baseEv + zeroEventIdx;
                    acc = acc + single(beamform_planewave_gpu(single(VB(ei).Data), ...
                        rxP, AP(zi).angle, xG, zG, fs, cc, doff, dT{1}));
                else
                    compound = complex(zeros(nZ,nX,'single'));
                    for a = 1:numAngles
                        rxPa = eP(AP(a).rxElements);
                        posEv = baseEv + AP(a).posIdx;
                        rfA = single(VB(posEv).Data);
                        if hasPulseInversion && ~isempty(AP(a).negIdx)
                            negEv = baseEv + AP(a).negIdx;
                            rfA = rfA + single(VB(negEv).Data);
                        end
                        compound = compound + single(beamform_planewave_gpu( ...
                            rfA, rxPa, AP(a).angle, xG, zG, fs, cc, doff, dT{a}));
                    end
                    acc = acc + compound / numAngles;
                end
            end
            clear VB; bgIQ = acc/nBF;
            save(bgCacheFile, 'bgIQ', '-v7.3');
            fprintf('    Background: %d frames, cached.\n', nBF);
        catch ME
            fprintf('    WARNING: Background failed: %s\n', ME.message);
        end
    end
end

totalFrames = 0;
for iB = 1:numel(blocks)
    cacheFile = fullfile(cacheDir, sprintf('block_%02d_IQ.mat', iB));
    
    if exist(cacheFile, 'file')
        fprintf('  Block %d/%d: cached (skipping beamform)\n', iB, numel(blocks));
        tmp = load(cacheFile, 'nFrames'); 
        totalFrames = totalFrames + tmp.nFrames;
        continue;
    end
    
    bt = tic;
    tail = blocks{iB}; if numel(tail)>35, tail=['...' tail(end-32:end)]; end
    fprintf('  Block %d/%d: %s\n', iB, numel(blocks), tail);
    
    [~,~,~,BC] = VsiVadaDataRead(dataFolder, blocks{iB}, ...
        1:fixed.eventsPerFrame, fixed.modeName);
    nEv = numel(BC.PulseSequences(1).Events);
    nCF = floor(nEv/fixed.eventsPerFrame);
    
    IQ_raw = zeros(nZ, nX, nCF, 'single');
    nCh = ceil(nCF/fixed.chunkSize);
    
    for iC = 1:nCh
        f1=(iC-1)*fixed.chunkSize+1; f2=min(iC*fixed.chunkSize,nCF); nF=f2-f1+1;
        e1=(f1-1)*fixed.eventsPerFrame+1; e2=f2*fixed.eventsPerFrame;
        try [VC,~,~,~]=VsiVadaDataRead(dataFolder,blocks{iB},e1:e2,fixed.modeName);
        catch, continue; end
        for f=1:nF
            baseEv = (f-1)*fixed.eventsPerFrame;
            if fixed.zeroOnly
                % Single zero-angle beamform
                ei = baseEv + zeroEventIdx;
                IQ_raw(:,:,f1+f-1) = single(beamform_planewave_gpu( ...
                    single(VC(ei).Data), rxP, AP(zi).angle, xG, zG, fs, cc, doff, dT{1}));
            else
                % Coherent compounding: beamform each angle, sum
                compound = complex(zeros(nZ, nX, 'single'));
                for a = 1:numAngles
                    rxPa = eP(AP(a).rxElements);
                    posEv = baseEv + AP(a).posIdx;
                    rfA = single(VC(posEv).Data);
                    if hasPulseInversion && ~isempty(AP(a).negIdx)
                        negEv = baseEv + AP(a).negIdx;
                        rfA = rfA + single(VC(negEv).Data);
                    end
                    compound = compound + single(beamform_planewave_gpu( ...
                        rfA, rxPa, AP(a).angle, xG, zG, fs, cc, doff, dT{a}));
                end
                IQ_raw(:,:,f1+f-1) = compound / numAngles;
            end
        end
        clear VC;
    end
    
    % Background subtraction
    if ~isempty(bgIQ)
        for f=1:nCF, IQ_raw(:,:,f)=IQ_raw(:,:,f)-bgIQ; end
    end
    
    nFrames = nCF; %#ok
    save(cacheFile, 'IQ_raw', 'nFrames', '-v7.3');
    totalFrames = totalFrames + nCF;
    fprintf('    %d frames, cached in %.1f min\n', nCF, toc(bt)/60);
    clear IQ_raw;

    if fixed.useGPU, try wait(gpuDevice); catch, end; end  % sync only — reset wipes delay tables
end

fprintf('  Total: %d frames across %d blocks\n', totalFrames, numel(blocks));
fprintf('  Beamforming: %.1f min\n', toc(bfTimer)/60);

%% ========================================================================
%  SVD FILTER + CACHE (per cutoff)
%  ========================================================================
fprintf('\n[3/5] SVD filtering at cutoffs %s...\n', mat2str(sweep.svdCutoffs));

for ci = 1:numel(sweep.svdCutoffs)
    cut = sweep.svdCutoffs(ci);
    svdCacheFile = fullfile(cacheDir, sprintf('svd_%d_filtered.mat', cut));
    
    if exist(svdCacheFile, 'file')
        fprintf('  SVD=%d: cached (skipping)\n', cut);
        continue;
    end
    
    fprintf('  SVD=%d: filtering all blocks...\n', cut);
    st = tic;
    
    locIdx = 0;
    gOff = 0; tF = 0;

    minThresh = min(sweep.detThresholds);
    detCfg.method = 'threshold';
    detCfg.threshold = minThresh;
    detCfg.minSep_mm = fixed.det.minSep_mm;
    detCfg.roiSize_px = fixed.det.roiSize_px;
    halfROI = floor(detCfg.roiSize_px / 2);

    % Preallocate locs in chunks to avoid O(N^2) reallocation
    chunkAlloc = 500000;
    allLocs = zeros(chunkAlloc, 5, 'single');  % [x, z, amplitude, globalFrameIdx, noiseStd]

    % Memory-efficient: load raw IQ one chunk at a time, SVD filter, detect,
    % then discard. Never hold a full block in RAM.
    svdChunkSize = 1500;  % Max frames per SVD pass (GPU memory safe)

    for iB = 1:numel(blocks)
        blockCache = fullfile(cacheDir, sprintf('block_%02d_IQ.mat', iB));

        % Get block frame count without loading full array
        info = whos('-file', blockCache, 'IQ_raw');
        nCF = info.size(3);
        nChunks = ceil(nCF / svdChunkSize);

        fprintf('    Block %d/%d: %d frames (%d SVD chunks)\n', iB, numel(blocks), nCF, nChunks);

        for iChk = 1:nChunks
            f1 = (iChk-1)*svdChunkSize + 1;
            f2 = min(iChk*svdChunkSize, nCF);
            nChkFrames = f2 - f1 + 1;

            % Load only this chunk from the cached block
            % matfile allows partial loading without reading full array
            mf = matfile(blockCache);
            IQ_chunk = mf.IQ_raw(:, :, f1:f2);

            fprintf('      Chunk %d/%d (frames %d-%d): SVD...', iChk, nChunks, f1, f2);
            tChk = tic;
            IQ_filt = svd_clutter_filter_gpu(gather(IQ_chunk), cut, [], fixed.useGPU);
            clear IQ_chunk;
            fprintf(' %.1fs. Detecting...', toc(tChk));

            % Detect frame-by-frame from this filtered chunk
            tDet = tic;
            nDetsChunk = 0;
            for f = 1:nChkFrames
                frame = abs(IQ_filt(:,:,f));
                noiseStd = std(frame(:));

                dets = detect_microbubbles(frame, detCfg, fixed.bf.dx, fixed.bf.dz);
                if isempty(dets), continue; end
                nDets = size(dets,1);
                nDetsChunk = nDetsChunk + nDets;

                % Ensure capacity
                if locIdx + nDets > size(allLocs,1)
                    allLocs(end+chunkAlloc, :) = 0;
                end

                for d = 1:nDets
                    r1 = max(1, dets(d,1)-halfROI); r2 = min(nZ, dets(d,1)+halfROI);
                    c1 = max(1, dets(d,2)-halfROI); c2 = min(nX, dets(d,2)+halfROI);
                    roi = frame(r1:r2, c1:c2);
                    [sR, sC] = intensity_weighted_centroid(roi);
                    locIdx = locIdx + 1;
                    allLocs(locIdx,:) = [xG(c1)+(sC-1)*fixed.bf.dx, ...
                        zG(r1)+(sR-1)*fixed.bf.dz, max(roi(:)), f1+f-1+gOff, noiseStd];
                end
            end
            clear IQ_filt;
            fprintf(' %d dets (%.1fs)\n', nDetsChunk, toc(tDet));
        end
        gOff = gOff + nCF;
    end

    allLocs = allLocs(1:locIdx, :);  % trim preallocated excess
    svdLocs = allLocs; %#ok
    save(svdCacheFile, 'svdLocs', 'gOff', '-v7.3');
    fprintf('    SVD=%d: %d locs from %d frames (%.1f min)\n', ...
        cut, size(allLocs,1), gOff, toc(st)/60);

    if fixed.useGPU, try wait(gpuDevice); catch, end; end  % sync only
end

%% ========================================================================
%  PARAMETER SWEEP
%  ========================================================================
fprintf('\n[4/5] LAT-ULM sweep: %d parameter combinations...\n', nCombos);

% Compute expected velocities for tube annotations
for i = 1:numel(tubes)
    area = pi * (tubes(i).ID_mm/2)^2;
    Q = tubes(i).rate_mL_min * 1000 / 60;  % mm^3/s
    tubes(i).expectedMean = Q / area;
    tubes(i).expectedPeak = 2 * tubes(i).expectedMean;
end

% Results storage
sweepResults = struct('svd',{}, 'threshold',{}, 'minTrackLen',{}, ...
    'maxGap',{}, 'maxDisp',{}, 'procNoise',{}, 'measNoise',{}, ...
    'nTracks',{}, 'nLocs',{}, 'medSpeed',{}, 'medLen',{}, 'folder',{});
comboIdx = 0;
comboTimer = tic;
sweepTimer = tic;

for ci = 1:numel(sweep.svdCutoffs)
    cut = sweep.svdCutoffs(ci);

    % Load cached locs for this SVD cutoff
    svdCacheFile = fullfile(cacheDir, sprintf('svd_%d_filtered.mat', cut));
    tmp = load(svdCacheFile, 'svdLocs', 'gOff');
    allLocsRaw = tmp.svdLocs;
    tF = tmp.gOff;

    for ti = 1:numel(sweep.detThresholds)
        thr = sweep.detThresholds(ti);

        % Filter localizations by threshold: keep only detections where
        % amplitude > thr * noiseStd (column 3 > thr * column 5)
        keepMask = allLocsRaw(:,3) > thr * allLocsRaw(:,5);
        locs = allLocsRaw(keepMask, 1:4);  % [x, z, amp, frameIdx]
        nLocs = size(locs, 1);

        % Skip if too many locs — noise-dominated, tracking would take forever
        maxLocsForTracking = 2e6;  % 2M locs max (~100 per frame avg for 66k frames is generous)
        if nLocs > maxLocsForTracking
            fprintf('  [SVD=%d, thr=%d] Skipping: %d locs > %.0fM limit (noise-dominated)\n', ...
                cut, thr, nLocs, maxLocsForTracking/1e6);
            % Fill sweep results with NaN for skipped combos
            for li_ = 1:numel(sweep.minTrackLengths)
                for gi_ = 1:numel(sweep.maxGapFrames)
                    for di_ = 1:numel(sweep.maxDisps)
                        for pi_ = 1:numel(sweep.processNoises)
                            for mi_ = 1:numel(sweep.measNoises)
                                comboIdx = comboIdx + 1;
                                sweepResults(comboIdx).svd = cut;
                                sweepResults(comboIdx).threshold = thr;
                                sweepResults(comboIdx).minTrackLen = sweep.minTrackLengths(li_);
                                sweepResults(comboIdx).maxGap = sweep.maxGapFrames(gi_);
                                sweepResults(comboIdx).maxDisp = sweep.maxDisps(di_);
                                sweepResults(comboIdx).procNoise = sweep.processNoises(pi_);
                                sweepResults(comboIdx).measNoise = sweep.measNoises(mi_);
                                sweepResults(comboIdx).nTracks = 0;
                                sweepResults(comboIdx).nLocs = nLocs;
                                sweepResults(comboIdx).medSpeed = NaN;
                                sweepResults(comboIdx).medLen = NaN;
                                sweepResults(comboIdx).folder = '';
                                sweepResults(comboIdx).skipped = true;
                            end
                        end
                    end
                end
            end
            continue;  % next threshold
        end

        for li = 1:numel(sweep.minTrackLengths)
            minLen = sweep.minTrackLengths(li);

            for gi = 1:numel(sweep.maxGapFrames)
                maxGap = sweep.maxGapFrames(gi);

                for di = 1:numel(sweep.maxDisps)
                    maxD = sweep.maxDisps(di);

                    for pi = 1:numel(sweep.processNoises)
                        pNoise = sweep.processNoises(pi);

                    for mi = 1:numel(sweep.measNoises)
                        mNoise = sweep.measNoises(mi);
                        comboIdx = comboIdx + 1;

            % Hierarchical output:
            %   LAT_ULM/svd2/thr5/len10_gap3_disp0.050_pn0.0005_mn0.080/
            comboLeaf = sprintf('len%d_gap%d_disp%.3f_pn%.4f_mn%.3f', ...
                minLen, maxGap, maxD, pNoise, mNoise);
            comboName = sprintf('svd%d/thr%d/%s', cut, thr, comboLeaf);
            comboDir = fullfile(outputBase, 'LAT_ULM', ...
                sprintf('svd%d', cut), sprintf('thr%d', thr), comboLeaf);

            % --- Skip if results already exist ---
            existingResult = fullfile(comboDir, 'LAT_ULM_results.mat');
            existingParams = fullfile(comboDir, 'params.txt');
            if exist(existingResult, 'file') || exist(existingParams, 'file')
                % Reload stats from params.txt if possible
                nTracks_loaded = 0; medSpeed_loaded = NaN; medLen_loaded = NaN;
                if exist(existingParams, 'file')
                    ptxt = fileread(existingParams);
                    tok = regexp(ptxt, 'Tracks:\s+(\d+)', 'tokens');
                    if ~isempty(tok), nTracks_loaded = str2double(tok{1}{1}); end
                    tok = regexp(ptxt, 'Median speed:\s+([\d.]+)', 'tokens');
                    if ~isempty(tok), medSpeed_loaded = str2double(tok{1}{1}); end
                    tok = regexp(ptxt, 'Median length:\s+([\d.]+)', 'tokens');
                    if ~isempty(tok), medLen_loaded = str2double(tok{1}{1}); end
                end
                sweepResults(comboIdx).svd = cut;
                sweepResults(comboIdx).threshold = thr;
                sweepResults(comboIdx).minTrackLen = minLen;
                sweepResults(comboIdx).maxGap = maxGap;
                sweepResults(comboIdx).maxDisp = maxD;
                sweepResults(comboIdx).procNoise = pNoise;
                sweepResults(comboIdx).measNoise = mNoise;
                sweepResults(comboIdx).nTracks = nTracks_loaded;
                sweepResults(comboIdx).nLocs = nLocs;
                sweepResults(comboIdx).medSpeed = medSpeed_loaded;
                sweepResults(comboIdx).medLen = medLen_loaded;
                sweepResults(comboIdx).folder = comboDir;
                sweepResults(comboIdx).skipped = false;
                if mod(comboIdx, 200) == 0
                    fprintf('  [%d/%d] %s — already done (skip)\n', comboIdx, nCombos, comboName);
                end
                continue;  % skip to next combo
            end

            if ~exist(comboDir, 'dir'), mkdir(comboDir); end

            if comboIdx > 1
                avgTime = toc(comboTimer) / (comboIdx - 1);
                etaMin = avgTime * (nCombos - comboIdx + 1) / 60;
                fprintf('  [%d/%d] %s (%d locs) [ETA: %.0f min]...', comboIdx, nCombos, comboName, nLocs, etaMin);
            else
                fprintf('  [%d/%d] %s (%d locs)...', comboIdx, nCombos, comboName, nLocs);
            end

            % Track
            trkP.maxDisp_mm     = maxD;
            trkP.maxGapFrames   = maxGap;
            trkP.minTrackLength = minLen;
            trkP.kalman.processNoise = pNoise;
            trkP.kalman.measNoise    = mNoise;

            tracks = track_microbubbles(locs, trkP, (0:tF-1)'/frameRate*1000);
            nTracks = numel(tracks);

            % Speed stats
            medSpeed = 0; medLen = 0;
            if nTracks > 0
                lens = cellfun(@(t) size(t,1), tracks);
                medLen = median(lens);
                speeds = cell(nTracks, 1);
                for iT = 1:nTracks
                    t = tracks{iT};
                    if size(t,1) > 1
                        dx=diff(t(:,1)); dz=diff(t(:,2)); df=diff(t(:,4));
                        dt=df/frameRate;
                        speeds{iT} = sqrt(dx.^2+dz.^2)./(dt+eps);
                    end
                end
                speeds = vertcat(speeds{:});
                if ~isempty(speeds), medSpeed = median(speeds); end
            end

            fprintf(' %d tracks (medSpd=%.1f mm/s, medLen=%.0f)\n', nTracks, medSpeed, medLen);

            % Store results
            sweepResults(comboIdx).svd = cut;
            sweepResults(comboIdx).threshold = thr;
            sweepResults(comboIdx).minTrackLen = minLen;
            sweepResults(comboIdx).maxGap = maxGap;
            sweepResults(comboIdx).maxDisp = maxD;
            sweepResults(comboIdx).procNoise = pNoise;
            sweepResults(comboIdx).measNoise = mNoise;
            sweepResults(comboIdx).nTracks = nTracks;
            sweepResults(comboIdx).nLocs = nLocs;
            sweepResults(comboIdx).medSpeed = medSpeed;
            sweepResults(comboIdx).medLen = medLen;
            sweepResults(comboIdx).folder = comboDir;

            % Write params.txt for easy reference
            fid = fopen(fullfile(comboDir, 'params.txt'), 'w');
            fprintf(fid, 'Dataset:        %s\n', label);
            fprintf(fid, 'SVD cutoff:     %d\n', cut);
            fprintf(fid, 'Threshold:      %d\n', thr);
            fprintf(fid, 'minTrackLength: %d\n', minLen);
            fprintf(fid, 'maxGapFrames:   %d\n', maxGap);
            fprintf(fid, 'maxDisp_mm:     %.3f\n', maxD);
            fprintf(fid, 'processNoise:   %.4f\n', pNoise);
            fprintf(fid, 'measNoise:      %.3f\n', mNoise);
            fprintf(fid, 'minSep_mm:      %.3f\n', fixed.det.minSep_mm);
            fprintf(fid, 'minSpeed_mm_s:  %.1f\n', fixed.postproc.minSpeed_mm_s);
            fprintf(fid, '--- Results ---\n');
            fprintf(fid, 'Tracks:         %d\n', nTracks);
            fprintf(fid, 'Localizations:  %d\n', nLocs);
            fprintf(fid, 'Median speed:   %.1f mm/s\n', medSpeed);
            fprintf(fid, 'Median length:  %.0f frames\n', medLen);
            fclose(fid);

            % Save results .mat (only for combos with tracks to save disk)
            if nTracks > 0
                R.localizations = locs;
                R.tracks = tracks;
                R.totalFrames = tF;
                R.config = fixed;
                R.config.svd.cutoffLow = cut;
                R.config.det.threshold = thr;
                R.config.track.minTrackLength = minLen;
                R.config.track.maxGapFrames = maxGap;
                R.config.track.maxDisp_mm = maxD;
                R.config.track.kalman.processNoise = pNoise;
                R.config.track.kalman.measNoise = mNoise;
                R.xGrid = xG; R.zGrid = zG; R.frameRate = frameRate;
                save(fullfile(comboDir, 'LAT_ULM_results.mat'), 'R', '-v7.3');
                clear R;
                % NOTE: Figures deferred to a separate figure step to avoid
                % CEF renderer crashes during overnight runs.
            else
                fprintf('    No tracks — skipping.\n');
            end

            clear tracks;

            % Checkpoint every 50 combos
            if mod(comboIdx, 50) == 0
                save(fullfile(outputBase, 'sweep_summary_partial.mat'), 'sweepResults', 'sweep', 'fixed');
                fprintf('  [Checkpoint] %d/%d combos saved\n', comboIdx, nCombos);
            end

                    end  % measNoise
                    end  % processNoise
                end  % maxDisp
            end  % maxGapFrames
        end  % minTrackLength
    end  % threshold

    % Intermediate save after each SVD cutoff
    if comboIdx > 0
        save(fullfile(outputBase, 'sweep_summary_partial.mat'), 'sweepResults', 'sweep', 'fixed');
        fprintf('  [Checkpoint] SVD=%d done — %d/%d combos, saved sweep_summary_partial.mat\n', cut, comboIdx, nCombos);
    end
end  % svdCutoff

elapsedTotal = toc(sweepTimer);
fprintf('\n  Sweep complete: %d combos in %.1f min\n', nCombos, elapsedTotal/60);

%% ========================================================================
%  SUMMARY COMPARISON
%  ========================================================================
fprintf('\n[Summary] Generating comparison...\n');

% Save sweep results
save(fullfile(outputBase, 'sweep_summary.mat'), 'sweepResults', 'sweep', 'fixed');

% Sort by track count (descending)
[~, sortIdx] = sort([sweepResults.nTracks], 'descend');

% Write full sweep_log.txt (human-readable, sorted by track count)
logFile = fullfile(outputBase, 'sweep_log.txt');
fid = fopen(logFile, 'w');
fprintf(fid, 'PARAMETER SWEEP RESULTS — %s\n', label);
fprintf(fid, 'Date: %s\n', datestr(now, 'yyyy-mm-dd HH:MM'));
fprintf(fid, 'Blocks: %d | Total frames: %d | Frame rate: %.0f Hz\n', numel(blocks), totalFrames, frameRate);
fprintf(fid, 'Total combos: %d\n\n', numel(sweepResults));
fprintf(fid, 'Sorted by track count (descending):\n');
fprintf(fid, '%-6s %-4s %-4s %-4s %-6s %-10s %-10s %7s %8s %8s  %s\n', ...
    'SVD', 'Thr', 'Len', 'Gap', 'Disp', 'ProcNoise', 'MeasNoise', 'Tracks', 'MedSpd', 'MedLen', 'Folder');
fprintf(fid, '%s\n', repmat('-', 1, 110));
for k = 1:numel(sweepResults)
    i = sortIdx(k);
    fprintf(fid, '%-6d %-4d %-4d %-4d %-6.3f %-10.4f %-10.3f %7d %8.1f %8.0f  %s\n', ...
        sweepResults(i).svd, sweepResults(i).threshold, sweepResults(i).minTrackLen, ...
        sweepResults(i).maxGap, sweepResults(i).maxDisp, sweepResults(i).procNoise, ...
        sweepResults(i).measNoise, sweepResults(i).nTracks, sweepResults(i).medSpeed, ...
        sweepResults(i).medLen, sweepResults(i).folder);
end
fclose(fid);
fprintf('  Saved: sweep_log.txt\n');

% Print top 30 to console
nShow = min(30, numel(sweepResults));
fprintf('\n  Top %d combos by track count:\n', nShow);
fprintf('  %-6s %-4s %-4s %-4s %-6s %-10s %-10s %7s %8s %8s\n', ...
    'SVD', 'Thr', 'Len', 'Gap', 'Disp', 'ProcNoise', 'MeasNoise', 'Tracks', 'MedSpd', 'MedLen');
fprintf('  %s\n', repmat('-', 1, 82));
for k = 1:nShow
    i = sortIdx(k);
    fprintf('  %-6d %-4d %-4d %-4d %-6.3f %-10.4f %-10.3f %7d %8.1f %8.0f\n', ...
        sweepResults(i).svd, sweepResults(i).threshold, sweepResults(i).minTrackLen, ...
        sweepResults(i).maxGap, sweepResults(i).maxDisp, sweepResults(i).procNoise, ...
        sweepResults(i).measNoise, sweepResults(i).nTracks, sweepResults(i).medSpeed, ...
        sweepResults(i).medLen);
end

% NOTE: Comparison figures deferred to a separate figure step
% Run that script after the sweep completes to generate all figures
% from saved .mat results — avoids CEF renderer crashes.

%% ========================================================================
%  [5/5] QUASAR SWEEP (per SVD cutoff × lambda)
%  ========================================================================
if ds.runQUASAR && exist('sushi_sparse_recovery','file') && exist('build_sushi_psf','file')
    fprintf('\n[5/5] QUASAR sweep...\n');
    qBlocks = ds.quasarBlocks;
    if isempty(qBlocks), qBlocks = blocks(1:min(2,numel(blocks))); end

    lambdas = defaults.quasar.lambda;
    srPx = defaults.quasar.srPixel_mm;
    ensSize = defaults.quasar.ensembleSize;

    for ci = 1:numel(sweep.svdCutoffs)
        cut = sweep.svdCutoffs(ci);
        for li = 1:numel(lambdas)
            lam = lambdas(li);
            qName = sprintf('svd%d/lambda%.2f', cut, lam);
            qDir = fullfile(outputBase, 'QUASAR', sprintf('svd%d', cut), ...
                sprintf('lambda%.2f', lam));
            if ~exist(qDir, 'dir'), mkdir(qDir); end

            % Check if already done
            if exist(fullfile(qDir, 'QUASAR_results.mat'), 'file')
                fprintf('  %s: cached (skipping)\n', qName);
                continue;
            end

            fprintf('  [QUASAR %d/%d] %s: processing %d blocks...\n', (ci-1)*numel(lambdas)+li, numel(sweep.svdCutoffs)*numel(lambdas), qName, numel(qBlocks));
            qt = tic;

            % Load and SVD-filter QUASAR blocks from cache
            qIQ = {};
            for qb = 1:numel(qBlocks)
                % Find block index in main block list
                bIdx = find(strcmp(blocks, qBlocks{qb}));
                if isempty(bIdx)
                    fprintf('    WARNING: QUASAR block not in main list, skipping\n');
                    continue;
                end
                blockCache = fullfile(cacheDir, sprintf('block_%02d_IQ.mat', bIdx));
                if ~exist(blockCache, 'file'), continue; end
                tmp = load(blockCache, 'IQ_raw');
                IQ_filt = svd_clutter_filter_gpu(tmp.IQ_raw, cut, [], fixed.useGPU);
                qIQ{end+1} = IQ_filt; %#ok
                clear tmp IQ_filt;
            end

            if isempty(qIQ)
                fprintf('    No valid QUASAR blocks — skipping\n');
                continue;
            end

            % Concatenate filtered IQ
            allQIQ = cat(3, qIQ{:}); clear qIQ;
            nQF = size(allQIQ, 3);

            % SR grid
            srX = xG(1):srPx:xG(end);
            srZ = zG(1):srPx:zG(end);

            % Build PSF — needs freq_MHz, c [m/s], pitch_mm, nRx, srPixel_mm, [nZ,nX], angle_deg
            nRxElements = numel(rxP);
            sNZ = numel(srZ); sNX = numel(srX);
            [psf, psfParams] = build_sushi_psf(fs, cc, pm, nRxElements, ...
                srPx, [sNZ, sNX], 0);

            % SUSHI + QUASAR opts
            fOpts.maxIter = 100;
            fOpts.nonNeg  = true;
            fOpts.useGPU  = fixed.useGPU;
            qOpts.maxIterCG     = 50;
            qOpts.supportThresh = 0;
            qOpts.useGPU        = fixed.useGPU;
            hasQuasarRefit = exist('quasar_refit', 'file') == 2;

            % Process ensembles (non-overlapping)
            nEns = floor(nQF / ensSize);
            sushiDensity  = zeros(sNZ, sNX, 'single');
            quasarDensity = zeros(sNZ, sNX, 'single');

            for e = 1:nEns
                idx = (e-1)*ensSize + (1:ensSize);
                ensIQ = allQIQ(:,:,idx);

                % Power image from ensemble, upscale to SR grid
                pw = mean(abs(ensIQ).^2, 3);
                clear ensIQ;
                pSR = max(imresize(pw, [sNZ, sNX], 'bicubic'), 0);
                nF = max(pSR(:)) + eps;

                % SUSHI sparse recovery (L1)
                try
                    [xf, ~] = sushi_sparse_recovery(pSR/nF, psf, lam, 'fista', fOpts);
                    sushiDensity = sushiDensity + xf * nF;

                    % QUASAR debiased refit on SUSHI support
                    if hasQuasarRefit
                        [xq, ~] = quasar_refit(pSR/nF, psf, xf, qOpts);
                        quasarDensity = quasarDensity + xq * nF;
                    end
                catch ME
                    fprintf('    Ensemble %d error: %s\n', e, ME.message);
                end
            end
            clear allQIQ;

            % Save QUASAR results
            R.sushiDensity  = sushiDensity;
            R.quasarDensity = quasarDensity;
            R.nEnsembles    = nEns;
            R.config = fixed;
            R.config.svd.cutoffLow = cut;
            R.config.quasar.lambda = lam;
            R.frameRate = frameRate;
            R.srX = srX; R.srZ = srZ;
            R.psfParams = psfParams;
            save(fullfile(qDir, 'QUASAR_results.mat'), 'R', '-v7.3');
            clear R;

            fprintf('    %s: %d ensembles in %.1f min\n', qName, nEns, toc(qt)/60);
            if fixed.useGPU, try wait(gpuDevice); catch, end; end  % sync only — preserve delay tables
        end
    end
else
    if ds.runQUASAR
        fprintf('\n[5/5] QUASAR: skipped (missing sushi_sparse_recovery or build_sushi_psf)\n');
    end
end

fprintf('\n=== Dataset %d/%d complete: %s (%.1f min) ===\n', ...
    iDS, nDS, ds.name, toc(dsTimer)/60);
if fixed.useGPU, try reset(gpuDevice); catch, end; end

end  % === END DATASET LOOP ===

fprintf('\n################################################################\n');
fprintf('  ALL DATASETS COMPLETE — %.1f hours total\n', toc(masterTimer)/3600);
fprintf('################################################################\n');

% Figure generation has been moved to GENERATE_SWEEP_FIGURES.m
% Run it after the sweep completes to generate three-panel and comparison figures.

function base = merge_struct(base, override)
    % Recursively merge override fields into base
    fns = fieldnames(override);
    for i = 1:numel(fns)
        if isstruct(override.(fns{i})) && isfield(base, fns{i}) && isstruct(base.(fns{i}))
            base.(fns{i}) = merge_struct(base.(fns{i}), override.(fns{i}));
        else
            base.(fns{i}) = override.(fns{i});
        end
    end
end
