%% PARAMETER_TUNER.m
% Fast parameter tuning tool for LAT-ULM pipeline.
%
% Beamforms a sample block ONCE (~2-3 min), caches to disk, then sweeps
% detection/tracking parameters near-instantly on the cached IQ stack.
% Total time: ~4-5 min first run, <1 min on cached re-runs.
%
% WORKFLOW:
%   1. Beamform + cache (or load cache)
%   2. SVD cutoff sweep → pick best cutoff
%   3. Detection threshold/minSep sweep → pick best detection params
%   4. Tracking maxDisp/minTrackLength sweep → pick best tracking params
%   5. Summary comparison figures
%   6. Export copy-pasteable config
%
% Author: Eli Wirth-Apley / Sun Lab, Northeastern University, March 2026

clearvars; close all; clc;

%% ========================================================================
%  SECTION 1: CONFIGURATION
%  ========================================================================

% --- Paths ---
tuner.dataFolder      = 'C:\path\to\VADA_data';
tuner.modeName        = '.vada';
tuner.vadaScriptsPath = 'C:\path\to\VevoF2\RFExampleScriptsF2\Matlab';
tuner.outputFolder    = fullfile(tuner.dataFolder, 'Tuner_Results');

% --- Input mode ---
% 'single' : process one block specified by tuner.baseFilename
% 'folder' : auto-discover all .vada blocks in tuner.dataFolder (with selection)
tuner.mode         = 'folder';  % 'single' or 'folder'
tuner.baseFilename = '';        % Only used in 'single' mode (no extension)

% --- Sample size ---
tuner.numFramesPerBlock = 500;  % Frames per block (200-2000). More = better stats, slower first run.

% --- GPU ---
tuner.useGPU = true;

% --- Angle mode (match pipeline) ---
tuner.zeroOnly = true;

% --- Overrides ---
tuner.sosOverride   = 1540;
tuner.pitchOverride = 0.300;

% --- Background ---
tuner.bgFile   = '';
tuner.bgFolder = '';

% --- Motion correction (applied before caching) ---
tuner.motionCorrection.enable    = false;
tuner.motionCorrection.method    = 'phase_corr';
tuner.motionCorrection.refType   = 'rolling';
tuner.motionCorrection.refWindow = 10;
tuner.motionCorrection.maxShift  = 5;

% --- Beamforming ([] = auto) ---
tuner.bf.xRange = [];
tuner.bf.zRange = [];
tuner.bf.dx     = [];
tuner.bf.dz     = [];

% --- Blanking ---
tuner.blankSteering = true;
tuner.blankMargin   = 1.5;
tuner.blankVoltage  = [];

% --- Parameter sweep ranges ---
tuner.svdCutoffs       = [1, 2, 5, 10, 20, 50];
tuner.detThresholds    = [3, 5, 8, 10, 15, 20];
tuner.detMinSeps_mm    = [0.05, 0.1, 0.2];
tuner.trackMaxDisps_mm    = [0.05, 0.1, 0.2, 0.5, 1.0];
tuner.trackMinLens        = [3, 5, 10, 20];
tuner.trackProcessNoises  = [0.0005, 0.001, 0.005, 0.01, 0.05];
tuner.trackMeasNoises     = [0.01, 0.03, 0.05, 0.08, 0.15];

% --- Defaults for non-swept parameters ---
tuner.det.roiSize_px = 7;
tuner.track.maxGapFrames = 3;

%% ========================================================================
%  SECTION 2: SETUP (metadata + delay tables)
%  ========================================================================
addpath(genpath(tuner.vadaScriptsPath));
pipelineDir = fileparts(mfilename('fullpath'));
if isempty(pipelineDir), pipelineDir = pwd; end
addpath(pipelineDir);
if ~exist(tuner.outputFolder, 'dir'), mkdir(tuner.outputFolder); end

if tuner.useGPU
    if gpuDeviceCount == 0
        warning('No GPU found.'); tuner.useGPU = false;
    else
        g = gpuDevice; reset(g);
        fprintf('GPU: %s (%.1f GB VRAM)\n', g.Name, g.TotalMemory/1e9);
    end
end

fprintf('\n=== PARAMETER TUNER ===\n');

%% File discovery
switch lower(tuner.mode)
    case 'single'
        if isempty(tuner.baseFilename)
            error('tuner.mode is ''single'' but tuner.baseFilename is empty.');
        end
        fileList    = {tuner.baseFilename};
        fileFolders = {tuner.dataFolder};

    case 'folder'
        fprintf('Scanning for %s files in:\n  %s\n', tuner.modeName, tuner.dataFolder);
        vadaFiles = dir(fullfile(tuner.dataFolder, '**', ['*' tuner.modeName]));
        if isempty(vadaFiles), error('No %s files found.', tuner.modeName); end

        nFound = numel(vadaFiles);
        fileList    = cell(nFound, 1);
        fileFolders = cell(nFound, 1);
        fileSizesGB = zeros(nFound, 1);
        for i = 1:nFound
            fullName = vadaFiles(i).name;
            fileList{i}    = fullName(1:end-numel(tuner.modeName));
            fileFolders{i} = vadaFiles(i).folder;
            fileSizesGB(i) = vadaFiles(i).bytes / 1e9;
        end
        [~, uIdx] = unique(fileList, 'stable');
        fileList = fileList(uIdx); fileFolders = fileFolders(uIdx);
        fileSizesGB = fileSizesGB(uIdx); nFound = numel(fileList);

        fprintf('\n  %-4s  %-50s  %8s\n', '#', 'Filename', 'Size(GB)');
        fprintf('  %s\n', repmat('-', 1, 66));
        for i = 1:nFound
            fprintf('  %-4d  %-50s  %8.1f\n', i, fileList{i}, fileSizesGB(i));
        end

        fprintf('\nSelect blocks (all / 1,3,5 / 1-3 / q): ');
        blockInput = input('', 's');
        blockInput = strtrim(blockInput);
        if strcmpi(blockInput, 'q'), fprintf('Cancelled.\n'); return; end

        if ~strcmpi(blockInput, 'all')
            selectedBlocks = [];
            parts = strsplit(blockInput, ',');
            for p = 1:numel(parts)
                tok = strtrim(parts{p});
                rangeParts = strsplit(tok, '-');
                if numel(rangeParts) == 2
                    r1 = str2double(rangeParts{1}); r2 = str2double(rangeParts{2});
                    if ~isnan(r1) && ~isnan(r2), selectedBlocks = [selectedBlocks, r1:r2]; end %#ok
                else
                    val = str2double(tok);
                    if ~isnan(val), selectedBlocks = [selectedBlocks, val]; end %#ok
                end
            end
            selectedBlocks = unique(selectedBlocks);
            selectedBlocks = selectedBlocks(selectedBlocks >= 1 & selectedBlocks <= nFound);
            if isempty(selectedBlocks), fprintf('No valid blocks.\n'); return; end
            fileList    = fileList(selectedBlocks);
            fileFolders = fileFolders(selectedBlocks);
        end

    otherwise
        error('tuner.mode must be ''single'' or ''folder''.');
end

numBlocks = numel(fileList);
fprintf('\nWill process %d block(s), %d frames each.\n', numBlocks, tuner.numFramesPerBlock);

% Auto-detect event structure from first block
numProbe = 30;
[VadaProbe, Param, TxrParam, ~] = VsiVadaDataRead(...
    fileFolders{1}, fileList{1}, 1:numProbe, tuner.modeName);

probeAngles = zeros(numProbe, 1);
probePolar  = zeros(numProbe, 1);
for ev = 1:numProbe
    if isfield(VadaProbe(ev).TxDelay, 'angle')
        probeAngles(ev) = VadaProbe(ev).TxDelay.angle;
    end
    if isfield(VadaProbe(ev).Waveform, 'Channel') && ...
            isfield(VadaProbe(ev).Waveform.Channel(1), 'invert')
        probePolar(ev) = VadaProbe(ev).Waveform.Channel(1).invert;
    end
end

probeSig = probeAngles * 10 + probePolar;
eventsPerFrame = 0;
for ev = 2:numProbe
    if probeSig(ev) == probeSig(1) && ev > 2
        patLen = ev - 1;
        if ev + patLen - 1 <= numProbe
            if all(probeSig(ev:ev+patLen-1) == probeSig(1:patLen))
                eventsPerFrame = patLen; break;
            end
        end
    end
end
if eventsPerFrame == 0, eventsPerFrame = 6; end

firstFrameAngles = probeAngles(1:eventsPerFrame);
uAngles = unique(firstFrameAngles, 'stable');
numAngles = numel(uAngles);
hasPI = all(arrayfun(@(a) sum(firstFrameAngles==a) >= 2, uAngles));

fprintf('  Detected: %d events/frame, %d angles, PI=%s\n', eventsPerFrame, numAngles, string(hasPI));

% Build anglePairs
VadaTest = VadaProbe(1:eventsPerFrame);
clear VadaProbe;

angleOrder = zeros(eventsPerFrame, 1);
polarityOrder = zeros(eventsPerFrame, 1);
for ev = 1:eventsPerFrame
    if isfield(VadaTest(ev).TxDelay, 'angle'), angleOrder(ev) = VadaTest(ev).TxDelay.angle; end
    if isfield(VadaTest(ev).Waveform, 'Channel') && isfield(VadaTest(ev).Waveform.Channel(1), 'invert')
        polarityOrder(ev) = VadaTest(ev).Waveform.Channel(1).invert;
    end
end

uniqueAngles = unique(angleOrder, 'stable');
anglePairs = struct('angle',{},'posIdx',{},'negIdx',{},'rxElements',{});
for a = 1:numel(uniqueAngles)
    anglePairs(a).angle = uniqueAngles(a);
    idxs = find(angleOrder == uniqueAngles(a));
    if numel(idxs) >= 2
        if polarityOrder(idxs(1)) == 0
            anglePairs(a).posIdx = idxs(1); anglePairs(a).negIdx = idxs(2);
        else
            anglePairs(a).posIdx = idxs(2); anglePairs(a).negIdx = idxs(1);
        end
    else
        anglePairs(a).posIdx = idxs(1); anglePairs(a).negIdx = [];
    end
    anglePairs(a).rxElements = VadaTest(anglePairs(a).posIdx).Elements;
end

% Probe parameters
rawPitch = TxrParam.ArrayPitch;
if rawPitch == 0 || isnan(rawPitch)
    if ~isempty(tuner.pitchOverride)
        pitch_mm = tuner.pitchOverride;
        fprintf('  Pitch: %.4f mm (override, ArrayPitch=0)\n', pitch_mm);
    else
        % Try XML
        xmlPath = fullfile(fileFolders{1}, [fileList{1} tuner.modeName '.xml']);
        if exist(xmlPath, 'file')
            xmlText = fileread(xmlPath);
            tok = regexp(xmlText, '<parameter\s+name="Element-Pitch"\s+value="([^"]*)"', 'tokens');
            if ~isempty(tok)
                pitch_mm = str2double(tok{1}{1});
                fprintf('  Pitch: %.4f mm (from XML, ArrayPitch=0)\n', pitch_mm);
            else
                error('ArrayPitch=0, pitchOverride=[], and no Element-Pitch in XML. Set tuner.pitchOverride.');
            end
        else
            error('ArrayPitch=0 and pitchOverride=[]. Set tuner.pitchOverride manually.');
        end
    end
elseif rawPitch < 10, pitch_mm = rawPitch;
else, pitch_mm = rawPitch / 1000; end

fs_MHz = Param.SampleFreq;
depthOffset_mm = Param.DepthOffset;
c_xml = Param.SoSMedia;
if ~isempty(tuner.sosOverride), c = tuner.sosOverride;
else
    c = c_xml;
    if c == 0, c = 1540; fprintf('  WARNING: SoS=0 in metadata, using 1540 m/s fallback\n'); end
    fprintf('  SoS: %.0f m/s (from metadata)\n', c);
end

txFreq_MHz = 6;
if isfield(VadaTest(1).Waveform, 'Channel') && ~isempty(VadaTest(1).Waveform.Channel)
    if isfield(VadaTest(1).Waveform.Channel(1), 'frequency')
        txFreq_MHz = VadaTest(1).Waveform.Channel(1).frequency;
    end
end
lambda_mm = c * 1e-3 / txFreq_MHz;

elemPos_mm = ((1:TxrParam.ArrayNumElements) - (TxrParam.ArrayNumElements+1)/2) * pitch_mm;
nSamplesPerEvent = size(VadaTest(1).Data, 1);
clear VadaTest;

% Grid auto-scaling
rxElemPos = elemPos_mm(anglePairs(1).rxElements);
rxSpan = max(rxElemPos) - min(rxElemPos);
maxDataDepth_mm = depthOffset_mm + (nSamplesPerEvent / (fs_MHz * 2)) * (c * 1e-3);

if isempty(tuner.bf.dx), tuner.bf.dx = round(lambda_mm / 5, 4); end
if isempty(tuner.bf.dz), tuner.bf.dz = round(lambda_mm / 10, 4); end
if isempty(tuner.bf.xRange)
    tuner.bf.xRange = [min(rxElemPos) - rxSpan*0.2, max(rxElemPos) + rxSpan*0.2];
end
if isempty(tuner.bf.zRange)
    tuner.bf.zRange = [max(depthOffset_mm, 0.5), maxDataDepth_mm];
end

xGrid = tuner.bf.xRange(1):tuner.bf.dx:tuner.bf.xRange(2);
zGrid = tuner.bf.zRange(1):tuner.bf.dz:tuner.bf.zRange(2);
nX = numel(xGrid); nZ = numel(zGrid);
dx = tuner.bf.dx; dz = tuner.bf.dz;

fprintf('  Grid: %d x %d (dx=%.4f, dz=%.4f mm)\n', nX, nZ, dx, dz);

% Identify 0-degree angle (needed for background + zeroOnly beamforming)
zeroAngleIdx = find([anglePairs.angle] == 0);
if isempty(zeroAngleIdx), [~, zeroAngleIdx] = min(abs([anglePairs.angle])); end
zeroEventIdx = anglePairs(zeroAngleIdx).posIdx;

% Delay tables
if tuner.zeroOnly
    rxPos_mm = elemPos_mm(anglePairs(zeroAngleIdx).rxElements);
    delayTables = {beamform_planewave_gpu([], rxPos_mm, anglePairs(zeroAngleIdx).angle, ...
        xGrid, zGrid, fs_MHz, c, depthOffset_mm, [])};
else
    delayTables = cell(numAngles, 1);
    for a = 1:numAngles
        rxPos_mm = elemPos_mm(anglePairs(a).rxElements);
        delayTables{a} = beamform_planewave_gpu([], rxPos_mm, anglePairs(a).angle, ...
            xGrid, zGrid, fs_MHz, c, depthOffset_mm, []);
    end
end

%% ========================================================================
%  SECTION 3: BEAMFORM + DISK CACHE (per-block, then concatenate)
%  ========================================================================

% Build a settings hash for cache invalidation. Any change to these
% parameters produces a different hash, forcing re-beamforming.
bgFileTag = tuner.bgFile; if isempty(bgFileTag), bgFileTag = 'none'; end
cacheSettings = sprintf('sos%.0f_pitch%.4f_dx%.4f_dz%.4f_x%.1f-%.1f_z%.1f-%.1f_zero%d_mc%d_nf%d_bg%s', ...
    c, pitch_mm, dx, dz, ...
    tuner.bf.xRange(1), tuner.bf.xRange(2), ...
    tuner.bf.zRange(1), tuner.bf.zRange(2), ...
    tuner.zeroOnly, tuner.motionCorrection.enable, tuner.numFramesPerBlock, bgFileTag);
cacheHash = dec2hex(mod(sum(double(cacheSettings) .* (1:numel(cacheSettings))), 2^32), 8);

fprintf('  Cache hash: %s (settings: %s)\n', cacheHash, cacheSettings);

% Compute background once (shared across blocks), with optional motion correction
bgMeanIQ = [];
if ~isempty(tuner.bgFile)
    bgDir = tuner.bgFolder; if isempty(bgDir), bgDir = tuner.dataFolder; end
    bgCacheFile = fullfile(tuner.outputFolder, sprintf('tuner_bg_%s_%s.mat', tuner.bgFile, cacheHash));

    if exist(bgCacheFile, 'file')
        fprintf('\n[Background] Loading cached background from:\n  %s\n', bgCacheFile);
        loaded = load(bgCacheFile, 'bgMeanIQ');
        bgMeanIQ = loaded.bgMeanIQ;
        fprintf('  Loaded background (motion-corrected=%s)\n', ...
            string(tuner.motionCorrection.enable));
    else
        fprintf('\n[Background] Computing from: %s\n', tuner.bgFile);
        try
            [VadaBg,~,~,BgCfg] = VsiVadaDataRead(bgDir, tuner.bgFile, ...
                1:eventsPerFrame, tuner.modeName);
            numBgEvents = numel(BgCfg.PulseSequences(1).Events);
            numBgFramesAvail = floor(numBgEvents / eventsPerFrame);
            numBgFrames = min(numBgFramesAvail, 500);
            clear VadaBg;

            fprintf('  Beamforming %d background frames...\n', numBgFrames);
            [VadaBg,~,~,~] = VsiVadaDataRead(bgDir, tuner.bgFile, ...
                1:numBgFrames*eventsPerFrame, tuner.modeName);

            % Beamform into full IQ stack (needed for motion correction)
            bgStack = zeros(nZ, nX, numBgFrames, 'single');
            for bf = 1:numBgFrames
                baseEvt = (bf-1) * eventsPerFrame;
                if tuner.zeroOnly
                    bgEvIdx = baseEvt + zeroEventIdx;
                    bgRf = single(VadaBg(bgEvIdx).Data);
                    rxPosBg = elemPos_mm(anglePairs(zeroAngleIdx).rxElements);
                    bgBf = beamform_planewave_gpu(bgRf, rxPosBg, ...
                        anglePairs(zeroAngleIdx).angle, ...
                        xGrid, zGrid, fs_MHz, c, depthOffset_mm, delayTables{1});
                else
                    bgBf = complex(zeros(nZ, nX, 'single'));
                    for a = 1:numAngles
                        rfPos = single(VadaBg(baseEvt + anglePairs(a).posIdx).Data);
                        if hasPI
                            rfNeg = single(VadaBg(baseEvt + anglePairs(a).negIdx).Data);
                            rfData = rfPos + rfNeg;
                        else
                            rfData = rfPos;
                        end
                        rxPosA = elemPos_mm(anglePairs(a).rxElements);
                        bfA = beamform_planewave_gpu(rfData, rxPosA, anglePairs(a).angle, ...
                            xGrid, zGrid, fs_MHz, c, depthOffset_mm, delayTables{a});
                        bgBf = bgBf + single(bfA);
                    end
                end
                bgStack(:,:,bf) = single(bgBf);
            end
            clear VadaBg;

            % Motion correction on background stack (in vivo: tissue moves)
            if tuner.motionCorrection.enable
                fprintf('  Motion-correcting background...\n');
                mcP.method = tuner.motionCorrection.method;
                mcP.refType = tuner.motionCorrection.refType;
                mcP.refWindow = tuner.motionCorrection.refWindow;
                mcP.maxShift = tuner.motionCorrection.maxShift;
                [bgShifts, bgMcDiag] = estimate_tissue_motion(bgStack, mcP);
                bgStack = apply_motion_correction(bgStack, bgShifts, tuner.useGPU);
                fprintf('  BG motion: max=%.2f px, mean=%.2f px\n', ...
                    bgMcDiag.maxDisp_px, bgMcDiag.meanDisp_px);
            end

            % Average the (motion-corrected) stack
            bgMeanIQ = mean(bgStack, 3);
            clear bgStack;

            % Cache
            save(bgCacheFile, 'bgMeanIQ', '-v7.3');
            fprintf('  Background cached: %d frames (motion-corrected=%s)\n', ...
                numBgFrames, string(tuner.motionCorrection.enable));
        catch ME
            fprintf('  WARNING: Background failed: %s\n', ME.message);
        end
    end
end

% Process each block: load cache or beamform
fprintf('\n[Beamform] Processing %d block(s), up to %d frames each...\n', ...
    numBlocks, tuner.numFramesPerBlock);

blockIQs = cell(numBlocks, 1);
blockTimestamps = cell(numBlocks, 1);
blockFrameCounts = zeros(numBlocks, 1);

for iBlock = 1:numBlocks
    cacheFile = fullfile(tuner.outputFolder, sprintf('tuner_%s_%s.mat', fileList{iBlock}, cacheHash));

    if exist(cacheFile, 'file')
        fprintf('  Block %d/%d [%s]: loading cache...', iBlock, numBlocks, fileList{iBlock});
        loaded = load(cacheFile, 'IQ_raw', 'timestamps');
        blockIQs{iBlock} = loaded.IQ_raw;
        blockTimestamps{iBlock} = loaded.timestamps;
        blockFrameCounts(iBlock) = size(loaded.IQ_raw, 3);
        fprintf(' %d frames\n', blockFrameCounts(iBlock));
    else
        fprintf('  Block %d/%d [%s]: beamforming...\n', iBlock, numBlocks, fileList{iBlock});
        tBF = tic;

        % Determine available frames
        [~,~,~,BlockConfig] = VsiVadaDataRead(fileFolders{iBlock}, fileList{iBlock}, ...
            1:eventsPerFrame, tuner.modeName);
        numTotalEvents = numel(BlockConfig.PulseSequences(1).Events);
        numAvailFrames = floor(numTotalEvents / eventsPerFrame);
        numFramesBlock = min(tuner.numFramesPerBlock, numAvailFrames);
        numEventsLoad = numFramesBlock * eventsPerFrame;

        fprintf('    Loading %d events (%d frames)...\n', numEventsLoad, numFramesBlock);
        [VadaData,~,~,~] = VsiVadaDataRead(fileFolders{iBlock}, fileList{iBlock}, ...
            1:numEventsLoad, tuner.modeName);

        IQ_block = zeros(nZ, nX, numFramesBlock, 'single');
        ts_block = zeros(numFramesBlock, 1);

        for iFrame = 1:numFramesBlock
            baseEvt = (iFrame-1) * eventsPerFrame;
            ts_block(iFrame) = VadaData(baseEvt+1).Timestamp;

            if tuner.zeroOnly
                evIdx = baseEvt + zeroEventIdx;
                rfData = single(VadaData(evIdx).Data);
                rxPosZ = elemPos_mm(anglePairs(zeroAngleIdx).rxElements);
                bfImg = beamform_planewave_gpu(rfData, rxPosZ, anglePairs(zeroAngleIdx).angle, ...
                    xGrid, zGrid, fs_MHz, c, depthOffset_mm, delayTables{1});
                IQ_block(:,:,iFrame) = single(bfImg);
            else
                compImg = complex(zeros(nZ, nX, 'single'));
                for a = 1:numAngles
                    rfPos = single(VadaData(baseEvt + anglePairs(a).posIdx).Data);
                    if hasPI
                        rfNeg = single(VadaData(baseEvt + anglePairs(a).negIdx).Data);
                        rfData = rfPos + rfNeg;
                    else
                        rfData = rfPos;
                    end
                    rxPosA = elemPos_mm(anglePairs(a).rxElements);
                    bfImg = beamform_planewave_gpu(rfData, rxPosA, anglePairs(a).angle, ...
                        xGrid, zGrid, fs_MHz, c, depthOffset_mm, delayTables{a});
                    compImg = compImg + single(bfImg);
                end
                IQ_block(:,:,iFrame) = compImg;
            end
        end
        clear VadaData;

        % Background subtraction
        if ~isempty(bgMeanIQ)
            for f = 1:numFramesBlock
                IQ_block(:,:,f) = IQ_block(:,:,f) - bgMeanIQ;
            end
        end

        % Motion correction
        if tuner.motionCorrection.enable
            mcP.method = tuner.motionCorrection.method;
            mcP.refType = tuner.motionCorrection.refType;
            mcP.refWindow = tuner.motionCorrection.refWindow;
            mcP.maxShift = tuner.motionCorrection.maxShift;
            [mcShifts, mcDiag] = estimate_tissue_motion(IQ_block, mcP);
            IQ_block = apply_motion_correction(IQ_block, mcShifts, tuner.useGPU);
            fprintf('    Motion: max=%.2f px, mean=%.2f px\n', mcDiag.maxDisp_px, mcDiag.meanDisp_px);
        end

        % Save per-block cache
        IQ_raw = IQ_block; timestamps = ts_block; %#ok (for save)
        save(cacheFile, 'IQ_raw', 'timestamps', 'xGrid', 'zGrid', '-v7.3');
        fprintf('    Cached in %.1f sec (%d frames)\n', toc(tBF), numFramesBlock);

        blockIQs{iBlock} = IQ_block;
        blockTimestamps{iBlock} = ts_block;
        blockFrameCounts(iBlock) = numFramesBlock;
    end
end

% Concatenate timestamps (keep blockIQs separate for per-block SVD)
timestamps = vertcat(blockTimestamps{:});
numFrames = sum(blockFrameCounts);
frameRate = 1000 / median(diff(timestamps));
fprintf('\n  Combined: %d blocks, %d total frames (%.1f sec of data), %.0f Hz\n', ...
    numBlocks, numFrames, numFrames/frameRate, frameRate);

%% ========================================================================
%  SECTION 4: SVD CUTOFF SWEEP (per-block, then concatenate)
%  ========================================================================
% SVD must run per-block to fit in GPU memory (same as the real pipeline).
% Each block's Casorati matrix [nPixels x nFramesBlock] is tractable.
fprintf('\n[SVD Sweep] Testing cutoffs: %s (per-block SVD)\n', mat2str(tuner.svdCutoffs));
tSVD = tic;

nCutoffs = numel(tuner.svdCutoffs);
svdStats = struct('cutoff', {}, 'snr', {}, 'meanEnvMax', {}, 'noiseStd', {}, ...
    'meanEnvImg', {}, 'mipImg', {}, 'tempStdImg', {});

% Process each cutoff one at a time — only keep 2D summary images
% for the figure, not the full 3D stack. Avoids OOM.
bestSNR = -Inf;
bestSVDIdx = 1;

for ci = 1:nCutoffs
    cut = tuner.svdCutoffs(ci);
    fprintf('    cutoff=%d...', cut);

    % SVD filter per-block, accumulate 2D statistics, then discard 3D data
    meanEnvAccum = zeros(nZ, nX, 'single');
    mipAccum     = zeros(nZ, nX, 'single');  % Max intensity projection
    stdAccumSum  = zeros(nZ, nX, 'double');   % For Welford online variance
    stdAccumSq   = zeros(nZ, nX, 'double');
    totalFr = 0;

    for iBlock = 1:numBlocks
        filtBlock = svd_clutter_filter_gpu( ...
            blockIQs{iBlock}, cut, [], tuner.useGPU);
        envBlock = abs(filtBlock);
        nFr = size(filtBlock, 3);

        meanEnvAccum = meanEnvAccum + sum(envBlock, 3);
        mipAccum = max(mipAccum, max(envBlock, [], 3));

        % Accumulate for temporal std (per-block, then combine)
        stdAccumSum = stdAccumSum + double(sum(envBlock, 3));
        stdAccumSq  = stdAccumSq  + double(sum(envBlock.^2, 3));
        totalFr = totalFr + nFr;
        clear filtBlock envBlock;
    end
    meanEnv = meanEnvAccum / totalFr;
    % Temporal std via E[X^2] - E[X]^2
    tempStd = single(sqrt(max(stdAccumSq/totalFr - (stdAccumSum/totalFr).^2, 0)));

    noiseReg = meanEnv(1:min(20, nZ), :);
    svdStats(ci).cutoff = cut;
    svdStats(ci).snr = 20*log10(max(meanEnv(:)) / (std(noiseReg(:)) + eps));
    svdStats(ci).meanEnvMax = max(meanEnv(:));
    svdStats(ci).noiseStd = std(noiseReg(:));
    svdStats(ci).meanEnvImg = meanEnv;
    svdStats(ci).mipImg = mipAccum;
    svdStats(ci).tempStdImg = tempStd;

    fprintf(' SNR=%.1f dB\n', svdStats(ci).snr);

    if svdStats(ci).snr > bestSNR
        bestSNR = svdStats(ci).snr;
        bestSVDIdx = ci;
    end
end

bestCutoff = tuner.svdCutoffs(bestSVDIdx);
fprintf('  >> Best SVD cutoff: %d (SNR=%.1f dB)\n', bestCutoff, bestSNR);

% Now rebuild ONLY the best cutoff's full stack for downstream sweeps
fprintf('  Rebuilding best cutoff stack (cut=%d)...\n', bestCutoff);
filteredBlocks = cell(numBlocks, 1);
for iBlock = 1:numBlocks
    filteredBlocks{iBlock} = svd_clutter_filter_gpu( ...
        blockIQs{iBlock}, bestCutoff, [], tuner.useGPU);
end
IQ_svd = cat(3, filteredBlocks{:});
clear filteredBlocks blockIQs;
fprintf('  SVD sweep done in %.1f sec.\n', toc(tSVD));

% --- Figure A: SVD comparison (3 rows: mean envelope, MIP, temporal std) ---
figure('Name', 'SVD Cutoff Sweep', 'Position', [50 50 1800 900]);
nCols = nCutoffs;

for ci = 1:nCutoffs
    markerStr = ''; if ci == bestSVDIdx, markerStr = ' <<'; end

    % Row 1: Mean envelope (tissue + static structure)
    subplot(3, nCols, ci);
    me = svdStats(ci).meanEnvImg;
    imagesc(xGrid, zGrid, 20*log10(me / max(me(:)) + eps));
    axis image; colormap(gca, gray); colorbar; caxis([-60 0]);
    title(sprintf('cut=%d%s', tuner.svdCutoffs(ci), markerStr), 'FontSize', 9);
    if ci == 1, ylabel({'Mean Envelope'; 'Ax [mm]'}); else, ylabel(''); end

    % Row 2: Max intensity projection (highlights brightest bubble events)
    subplot(3, nCols, nCols + ci);
    mip = svdStats(ci).mipImg;
    imagesc(xGrid, zGrid, 20*log10(mip / max(mip(:)) + eps));
    axis image; colormap(gca, hot); colorbar; caxis([-40 0]);
    if ci == 1, ylabel({'MIP'; 'Ax [mm]'}); else, ylabel(''); end

    % Row 3: Temporal std (highlights moving scatterers = bubbles)
    subplot(3, nCols, 2*nCols + ci);
    ts = svdStats(ci).tempStdImg;
    imagesc(xGrid, zGrid, 20*log10(ts / max(ts(:)) + eps));
    axis image; colormap(gca, hot); colorbar; caxis([-40 0]);
    xlabel('Lat [mm]');
    if ci == 1, ylabel({'Temporal Std'; 'Ax [mm]'}); else, ylabel(''); end
end
sgtitle(sprintf('SVD Cutoff Sweep (%d frames, %d blocks) — Row 1: Mean Env, Row 2: MIP, Row 3: Temporal Std', ...
    numFrames, numBlocks));
saveas(gcf, fullfile(tuner.outputFolder, 'tuner_svd_sweep.png'));

% Print SVD stats
fprintf('\n  %-8s %-10s %-12s %-12s\n', 'Cutoff', 'SNR(dB)', 'MaxSignal', 'NoiseStd');
fprintf('  %s\n', repmat('-', 1, 42));
for ci = 1:nCutoffs
    marker = ''; if ci == bestSVDIdx, marker = ' <<'; end
    fprintf('  %-8d %-10.1f %-12.1f %-12.1f%s\n', svdStats(ci).cutoff, ...
        svdStats(ci).snr, svdStats(ci).meanEnvMax, svdStats(ci).noiseStd, marker);
end

%% ========================================================================
%  SECTION 5: DETECTION SWEEP
%  ========================================================================
fprintf('\n[Detection Sweep] Thresholds: %s, MinSeps: %s mm\n', ...
    mat2str(tuner.detThresholds), mat2str(tuner.detMinSeps_mm));
tDet = tic;

nThresh = numel(tuner.detThresholds);
nMinSep = numel(tuner.detMinSeps_mm);
detResults = zeros(nThresh, nMinSep);  % mean locs/frame

% --- FAST PASS: count detections only (no sub-pixel localization) ---
% Pre-compute all envelopes once
fprintf('  Pre-computing envelopes...\n');
envFrames = abs(IQ_svd);

detParams.method = 'threshold';
detParams.roiSize_px = tuner.det.roiSize_px;
halfROI = floor(detParams.roiSize_px / 2);

for si = 1:nMinSep
    detParams.minSep_mm = tuner.detMinSeps_mm(si);
    for ti = 1:nThresh
        detParams.threshold = tuner.detThresholds(ti);
        totalLocs = 0;
        for f = 1:numFrames
            dets = detect_microbubbles(envFrames(:,:,f), detParams, dx, dz);
            if ~isempty(dets), totalLocs = totalLocs + size(dets, 1); end
        end
        detResults(ti, si) = totalLocs / numFrames;
    end
    fprintf('    minSep=%.3f mm done\n', tuner.detMinSeps_mm(si));
end
fprintf('  Detection count sweep done in %.1f sec.\n', toc(tDet));

% --- FULL PASS: sub-pixel localization for best combo only ---
% Auto-select: closest to target locs/frame with middle minSep
midS = ceil(nMinSep/2);
[~, bestDetIdx] = min(abs(detResults(:, midS) - 3));
bestThreshold = tuner.detThresholds(bestDetIdx);
bestMinSep = tuner.detMinSeps_mm(midS);

fprintf('  Running sub-pixel localization for best combo (thresh=%d, minSep=%.3f)...\n', ...
    bestThreshold, bestMinSep);
tLoc = tic;

detParams.threshold = bestThreshold;
detParams.minSep_mm = bestMinSep;
bestLocs = zeros(0, 4);

for f = 1:numFrames
    frame = envFrames(:,:,f);
    dets = detect_microbubbles(frame, detParams, dx, dz);
    if isempty(dets), continue; end
    for d = 1:size(dets, 1)
        r1 = max(1, dets(d,1)-halfROI); r2 = min(nZ, dets(d,1)+halfROI);
        c1 = max(1, dets(d,2)-halfROI); c2 = min(nX, dets(d,2)+halfROI);
        roi = frame(r1:r2, c1:c2);
        [sR, sC] = intensity_weighted_centroid(roi);
        bestLocs(end+1,:) = [xGrid(c1)+(sC-1)*dx, zGrid(r1)+(sR-1)*dz, max(roi(:)), f]; %#ok
    end
end
fprintf('  Localization done in %.1f sec (%d locs).\n', toc(tLoc), size(bestLocs,1));
clear envFrames;

% --- Figure B: Threshold tuning curves ---
figure('Name', 'Detection Sweep', 'Position', [50 50 1600 500]);

subplot(1,3,1);
colors = lines(nMinSep);
hold on;
for si = 1:nMinSep
    semilogy(tuner.detThresholds, detResults(:,si), '.-', 'LineWidth', 2, ...
        'MarkerSize', 12, 'Color', colors(si,:), ...
        'DisplayName', sprintf('minSep=%.2f mm', tuner.detMinSeps_mm(si)));
end
yline(3, 'g--', 'Target: 1-5/frame');
yline(1, 'g--'); yline(5, 'g--');
xlabel('Threshold (noise std)'); ylabel('Localizations / Frame');
title('Threshold Tuning'); legend('Location', 'best'); grid on;

% --- Figure C: Sample frame with detections ---
subplot(1,3,2);
sampleFrame = round(numFrames / 2);
frame = abs(IQ_svd(:,:,sampleFrame));
imagesc(xGrid, zGrid, frame); axis image; colormap(gca, hot); colorbar; hold on;

% Show detections from best combo on sample frame
fLocs = bestLocs(bestLocs(:,4) == sampleFrame, :);
if ~isempty(fLocs)
    plot(fLocs(:,1), fLocs(:,2), 'go', 'MarkerSize', 8, 'LineWidth', 2);
end
title(sprintf('Frame %d (thresh=%d, minSep=%.3f)', sampleFrame, bestThreshold, bestMinSep));
xlabel('Lat [mm]'); ylabel('Ax [mm]');

% --- Figure D: Quick density for best params ---
subplot(1,3,3);
if ~isempty(bestLocs)
    quickDens = zeros(nZ, nX);
    for d = 1:size(bestLocs, 1)
        [~, xi] = min(abs(xGrid - bestLocs(d,1)));
        [~, zi] = min(abs(zGrid - bestLocs(d,2)));
        quickDens(zi, xi) = quickDens(zi, xi) + 1;
    end
    imagesc(xGrid, zGrid, log10(quickDens + 1));
    axis image; colormap(gca, hot); colorbar;
    title(sprintf('Density (%.1f locs/frame)', detResults(bestDetIdx, midS)));
else
    text(0.5, 0.5, 'No detections', 'HorizontalAlignment', 'center'); axis off;
end
xlabel('Lat [mm]'); ylabel('Ax [mm]');

sgtitle(sprintf('Detection Sweep (SVD cut=%d, %d frames)', bestCutoff, numFrames));
saveas(gcf, fullfile(tuner.outputFolder, 'tuner_detection_sweep.png'));

fprintf('\n  >> Best detection: threshold=%d, minSep=%.3f mm (%.1f locs/frame)\n', ...
    bestThreshold, bestMinSep, detResults(bestDetIdx, midS));

%% ========================================================================
%  SECTION 6a: TRACKING SWEEP — maxDisp x minTrackLength
%  ========================================================================
% First sweep geometry params (maxDisp, minLen) with mid-range Kalman noise.
% Then sweep Kalman params using the best geometry.
midPN = tuner.trackProcessNoises(ceil(numel(tuner.trackProcessNoises)/2));
midMN = tuner.trackMeasNoises(ceil(numel(tuner.trackMeasNoises)/2));

fprintf('\n[Tracking Sweep 6a] MaxDisps: %s mm, MinLens: %s\n', ...
    mat2str(tuner.trackMaxDisps_mm), mat2str(tuner.trackMinLens));
fprintf('  Using Kalman: processNoise=%.4f, measNoise=%.3f (mid-range)\n', midPN, midMN);
tTrack = tic;

nDisp = numel(tuner.trackMaxDisps_mm);
nLen  = numel(tuner.trackMinLens);
trackCounts  = zeros(nDisp, nLen);
trackMedLen  = zeros(nDisp, nLen);
trackMedSpd  = zeros(nDisp, nLen);
trackResults = cell(nDisp, nLen);

for di = 1:nDisp
    for li = 1:nLen
        trkParams.maxDisp_mm     = tuner.trackMaxDisps_mm(di);
        trkParams.maxGapFrames   = tuner.track.maxGapFrames;
        trkParams.minTrackLength = tuner.trackMinLens(li);
        trkParams.kalman.processNoise = midPN;
        trkParams.kalman.measNoise    = midMN;

        tracks = track_microbubbles(bestLocs, trkParams, timestamps);
        trackResults{di, li} = tracks;
        trackCounts(di, li) = numel(tracks);

        if ~isempty(tracks)
            lens = cellfun(@(t) size(t,1), tracks);
            trackMedLen(di, li) = median(lens);

            speeds = [];
            for iT = 1:numel(tracks)
                t = tracks{iT};
                if size(t,1) > 1
                    speeds = [speeds; sqrt(diff(t(:,1)).^2 + diff(t(:,2)).^2) * frameRate]; %#ok
                end
            end
            if ~isempty(speeds), trackMedSpd(di, li) = median(speeds); end
        end
    end
end
fprintf('  Geometry sweep done in %.1f sec.\n', toc(tTrack));

% --- Figure E: Track count heatmap ---
figure('Name', 'Tracking Sweep — Geometry', 'Position', [50 50 1600 800]);

subplot(2,3,1);
imagesc(tuner.trackMinLens, tuner.trackMaxDisps_mm, trackCounts);
colorbar; xlabel('minTrackLength'); ylabel('maxDisp [mm]');
title('Track Count'); set(gca, 'YDir', 'normal');
for di = 1:nDisp
    for li = 1:nLen
        text(tuner.trackMinLens(li), tuner.trackMaxDisps_mm(di), ...
            sprintf('%d', trackCounts(di,li)), 'HorizontalAlignment', 'center', ...
            'Color', 'w', 'FontWeight', 'bold', 'FontSize', 8);
    end
end

subplot(2,3,2);
imagesc(tuner.trackMinLens, tuner.trackMaxDisps_mm, trackMedLen);
colorbar; xlabel('minTrackLength'); ylabel('maxDisp [mm]');
title('Median Track Length'); set(gca, 'YDir', 'normal');

subplot(2,3,3);
imagesc(tuner.trackMinLens, tuner.trackMaxDisps_mm, trackMedSpd);
colorbar; xlabel('minTrackLength'); ylabel('maxDisp [mm]');
title('Median Speed [mm/s]'); set(gca, 'YDir', 'normal');

% Track plots for 3 combos
midL = ceil(nLen/2);
dispIdxs = [1, ceil(nDisp/2), nDisp];
for pi = 1:3
    subplot(2,3,3+pi);
    di = dispIdxs(pi);
    tracks = trackResults{di, midL};
    hold on;
    if ~isempty(tracks)
        nPlot = min(numel(tracks), 200);
        cmap = jet(nPlot);
        plotIdx = randperm(numel(tracks), nPlot);
        for ii = 1:nPlot
            t = tracks{plotIdx(ii)};
            plot(t(:,1), t(:,2), '-', 'Color', [cmap(ii,:) 0.5], 'LineWidth', 0.5);
        end
    end
    set(gca, 'YDir', 'reverse'); axis equal tight;
    xlabel('Lat [mm]'); ylabel('Ax [mm]');
    title(sprintf('maxDisp=%.2f, minLen=%d\n%d tracks', ...
        tuner.trackMaxDisps_mm(di), tuner.trackMinLens(midL), trackCounts(di, midL)), 'FontSize', 9);
end

sgtitle(sprintf('Tracking Geometry Sweep (thresh=%d, %d locs)', bestThreshold, size(bestLocs,1)));
saveas(gcf, fullfile(tuner.outputFolder, 'tuner_tracking_geometry.png'));

% Auto-select best geometry: most tracks with median length > 10
validMask = trackMedLen >= 10;
if any(validMask(:))
    maskedCounts = trackCounts;
    maskedCounts(~validMask) = 0;
    [~, bestIdx] = max(maskedCounts(:));
    [bestDI, bestLI] = ind2sub([nDisp, nLen], bestIdx);
else
    [~, bestIdx] = max(trackCounts(:));
    [bestDI, bestLI] = ind2sub([nDisp, nLen], bestIdx);
end
bestMaxDisp = tuner.trackMaxDisps_mm(bestDI);
bestMinLen  = tuner.trackMinLens(bestLI);
fprintf('\n  >> Best geometry: maxDisp=%.3f mm, minLen=%d (%d tracks, medLen=%.0f)\n', ...
    bestMaxDisp, bestMinLen, trackCounts(bestDI, bestLI), trackMedLen(bestDI, bestLI));

%% ========================================================================
%  SECTION 6b: TRACKING SWEEP — Kalman processNoise x measNoise
%  ========================================================================
fprintf('\n[Tracking Sweep 6b] ProcessNoises: %s, MeasNoises: %s\n', ...
    mat2str(tuner.trackProcessNoises), mat2str(tuner.trackMeasNoises));
fprintf('  Using best geometry: maxDisp=%.3f, minLen=%d\n', bestMaxDisp, bestMinLen);
tKalman = tic;

nPN = numel(tuner.trackProcessNoises);
nMN = numel(tuner.trackMeasNoises);
kalmanCounts  = zeros(nPN, nMN);
kalmanMedLen  = zeros(nPN, nMN);
kalmanMedSpd  = zeros(nPN, nMN);
kalmanResults = cell(nPN, nMN);

for pi = 1:nPN
    for mi = 1:nMN
        trkParams.maxDisp_mm     = bestMaxDisp;
        trkParams.maxGapFrames   = tuner.track.maxGapFrames;
        trkParams.minTrackLength = bestMinLen;
        trkParams.kalman.processNoise = tuner.trackProcessNoises(pi);
        trkParams.kalman.measNoise    = tuner.trackMeasNoises(mi);

        tracks = track_microbubbles(bestLocs, trkParams, timestamps);
        kalmanResults{pi, mi} = tracks;
        kalmanCounts(pi, mi) = numel(tracks);

        if ~isempty(tracks)
            lens = cellfun(@(t) size(t,1), tracks);
            kalmanMedLen(pi, mi) = median(lens);

            speeds = [];
            for iT = 1:numel(tracks)
                t = tracks{iT};
                if size(t,1) > 1
                    speeds = [speeds; sqrt(diff(t(:,1)).^2 + diff(t(:,2)).^2) * frameRate]; %#ok
                end
            end
            if ~isempty(speeds), kalmanMedSpd(pi, mi) = median(speeds); end
        end
    end
end
fprintf('  Kalman sweep done in %.1f sec.\n', toc(tKalman));

% --- Figure: Kalman parameter heatmaps ---
figure('Name', 'Tracking Sweep — Kalman', 'Position', [50 50 1600 800]);

subplot(2,3,1);
imagesc(tuner.trackMeasNoises, tuner.trackProcessNoises, kalmanCounts);
colorbar; xlabel('measNoise'); ylabel('processNoise');
title('Track Count'); set(gca, 'YDir', 'normal', 'XScale', 'log', 'YScale', 'log');
for pi = 1:nPN
    for mi = 1:nMN
        text(tuner.trackMeasNoises(mi), tuner.trackProcessNoises(pi), ...
            sprintf('%d', kalmanCounts(pi,mi)), 'HorizontalAlignment', 'center', ...
            'Color', 'w', 'FontWeight', 'bold', 'FontSize', 8);
    end
end

subplot(2,3,2);
imagesc(tuner.trackMeasNoises, tuner.trackProcessNoises, kalmanMedLen);
colorbar; xlabel('measNoise'); ylabel('processNoise');
title('Median Track Length'); set(gca, 'YDir', 'normal');

subplot(2,3,3);
imagesc(tuner.trackMeasNoises, tuner.trackProcessNoises, kalmanMedSpd);
colorbar; xlabel('measNoise'); ylabel('processNoise');
title('Median Speed [mm/s]'); set(gca, 'YDir', 'normal');

% Track plots for 3 Kalman combos: low/mid/high processNoise at mid measNoise
midMI = ceil(nMN/2);
pnIdxs = [1, ceil(nPN/2), nPN];
for pi = 1:3
    subplot(2,3,3+pi);
    pIdx = pnIdxs(pi);
    tracks = kalmanResults{pIdx, midMI};
    hold on;
    if ~isempty(tracks)
        nPlot = min(numel(tracks), 200);
        cmap = jet(nPlot);
        plotIdx = randperm(numel(tracks), nPlot);
        for ii = 1:nPlot
            t = tracks{plotIdx(ii)};
            plot(t(:,1), t(:,2), '-', 'Color', [cmap(ii,:) 0.5], 'LineWidth', 0.5);
        end
    end
    set(gca, 'YDir', 'reverse'); axis equal tight;
    xlabel('Lat [mm]'); ylabel('Ax [mm]');
    title(sprintf('pNoise=%.4f, mNoise=%.3f\n%d tracks', ...
        tuner.trackProcessNoises(pIdx), tuner.trackMeasNoises(midMI), ...
        kalmanCounts(pIdx, midMI)), 'FontSize', 9);
end

sgtitle(sprintf('Kalman Noise Sweep (maxDisp=%.3f, minLen=%d)', bestMaxDisp, bestMinLen));
saveas(gcf, fullfile(tuner.outputFolder, 'tuner_tracking_kalman.png'));

% Auto-select best Kalman: most tracks with median length > 10
validKalman = kalmanMedLen >= 10;
if any(validKalman(:))
    maskedK = kalmanCounts;
    maskedK(~validKalman) = 0;
    [~, bestKIdx] = max(maskedK(:));
    [bestPI, bestMI] = ind2sub([nPN, nMN], bestKIdx);
else
    [~, bestKIdx] = max(kalmanCounts(:));
    [bestPI, bestMI] = ind2sub([nPN, nMN], bestKIdx);
end
bestProcessNoise = tuner.trackProcessNoises(bestPI);
bestMeasNoise    = tuner.trackMeasNoises(bestMI);
fprintf('\n  >> Best Kalman: processNoise=%.4f, measNoise=%.3f (%d tracks, medLen=%.0f)\n', ...
    bestProcessNoise, bestMeasNoise, kalmanCounts(bestPI, bestMI), kalmanMedLen(bestPI, bestMI));

% Re-run tracking with all best params for downstream use
trkParamsBest.maxDisp_mm     = bestMaxDisp;
trkParamsBest.maxGapFrames   = tuner.track.maxGapFrames;
trkParamsBest.minTrackLength = bestMinLen;
trkParamsBest.kalman.processNoise = bestProcessNoise;
trkParamsBest.kalman.measNoise    = bestMeasNoise;

%% ========================================================================
%  SECTION 7: SUMMARY COMPARISON
%  ========================================================================
fprintf('\n[Summary] Generating comparison figure...\n');

bestTracks = track_microbubbles(bestLocs, trkParamsBest, timestamps);

figure('Name', 'Tuner Summary', 'Position', [50 50 1600 700]);

% Row 1: SVD frame, detections, density
subplot(2,3,1);
meanEnv = mean(abs(IQ_svd), 3);
imagesc(xGrid, zGrid, 20*log10(meanEnv / max(meanEnv(:)) + eps));
axis image; colormap(gca, gray); colorbar; caxis([-60 0]);
title(sprintf('SVD filtered (cut=%d)', bestCutoff));
xlabel('Lat [mm]'); ylabel('Ax [mm]');

subplot(2,3,2);
frame = abs(IQ_svd(:,:,sampleFrame));
imagesc(xGrid, zGrid, frame); axis image; colormap(gca, hot); colorbar; hold on;
fLocs = bestLocs(bestLocs(:,4) == sampleFrame, :);
if ~isempty(fLocs)
    plot(fLocs(:,1), fLocs(:,2), 'go', 'MarkerSize', 10, 'LineWidth', 2);
end
title(sprintf('Detections (thresh=%d)', bestThreshold));
xlabel('Lat [mm]'); ylabel('Ax [mm]');

subplot(2,3,3);
if ~isempty(bestLocs)
    quickDens = zeros(nZ, nX);
    for d = 1:size(bestLocs, 1)
        [~, xi] = min(abs(xGrid - bestLocs(d,1)));
        [~, zi] = min(abs(zGrid - bestLocs(d,2)));
        quickDens(zi, xi) = quickDens(zi, xi) + 1;
    end
    imagesc(xGrid, zGrid, log10(quickDens + 1));
    axis image; colormap(gca, hot); colorbar;
end
title(sprintf('Density (%d locs)', size(bestLocs,1)));
xlabel('Lat [mm]'); ylabel('Ax [mm]');

% Row 2: tracks, speed dist, detection rate
subplot(2,3,4);
hold on;
if ~isempty(bestTracks)
    nPlot = min(numel(bestTracks), 500);
    cmap = jet(nPlot);
    plotIdx = randperm(numel(bestTracks), nPlot);
    for ii = 1:nPlot
        t = bestTracks{plotIdx(ii)};
        plot(t(:,1), t(:,2), '-', 'Color', [cmap(ii,:) 0.5], 'LineWidth', 0.5);
    end
end
set(gca, 'YDir', 'reverse'); axis equal tight;
xlabel('Lat [mm]'); ylabel('Ax [mm]');
title(sprintf('Tracks (%d/%d shown)', min(numel(bestTracks),500), numel(bestTracks)));

subplot(2,3,5);
if ~isempty(bestTracks)
    speeds = [];
    for iT = 1:numel(bestTracks)
        t = bestTracks{iT};
        if size(t,1) > 1
            speeds = [speeds; sqrt(diff(t(:,1)).^2 + diff(t(:,2)).^2) * frameRate]; %#ok
        end
    end
    if ~isempty(speeds)
        histogram(speeds, 50, 'FaceColor', [0.3 0.5 0.8]);
        xlabel('Speed [mm/s]'); ylabel('Count');
        title(sprintf('Speed (median=%.1f mm/s)', median(speeds)));
    end
end

subplot(2,3,6);
if ~isempty(bestLocs)
    lpf = accumarray(bestLocs(:,4), 1, [numFrames 1]);
    plot(lpf, 'Color', [0.2 0.6 0.3]); hold on;
    plot(movmean(lpf, 50), 'r-', 'LineWidth', 2);
    xlabel('Frame'); ylabel('Detections');
    title(sprintf('Detection Rate (mean=%.1f/frame)', mean(lpf)));
end

sgtitle(sprintf('TUNER SUMMARY: SVD=%d, thresh=%d, minSep=%.2f, maxDisp=%.3f, minLen=%d, pN=%.4f, mN=%.3f', ...
    bestCutoff, bestThreshold, bestMinSep, bestMaxDisp, bestMinLen, bestProcessNoise, bestMeasNoise));
saveas(gcf, fullfile(tuner.outputFolder, 'tuner_summary.png'));

%% ========================================================================
%  SECTION 8: EXPORT RECOMMENDED CONFIG
%  ========================================================================
fprintf('\n');
fprintf('=========================================================================\n');
fprintf('  RECOMMENDED CONFIGURATION (copy-paste into LAT_ULM_pipeline.m)\n');
fprintf('=========================================================================\n');
fprintf('\n');
fprintf('config.svd.cutoffLow  = %d;\n', bestCutoff);
fprintf('config.svd.cutoffHigh = [];\n');
fprintf('config.det.threshold  = %d;\n', bestThreshold);
fprintf('config.det.minSep_mm  = %.3f;\n', bestMinSep);
fprintf('config.det.roiSize_px = %d;\n', tuner.det.roiSize_px);
fprintf('config.track.maxDisp_mm     = %.3f;\n', bestMaxDisp);
fprintf('config.track.maxGapFrames   = %d;\n', tuner.track.maxGapFrames);
fprintf('config.track.minTrackLength = %d;\n', bestMinLen);
fprintf('config.track.kalman.processNoise = %.4f;\n', bestProcessNoise);
fprintf('config.track.kalman.measNoise    = %.3f;\n', bestMeasNoise);
fprintf('\n');
fprintf('=========================================================================\n');

% Save tuned config
tunedConfig.svd.cutoffLow  = bestCutoff;
tunedConfig.svd.cutoffHigh = [];
tunedConfig.det.method     = 'threshold';
tunedConfig.det.threshold  = bestThreshold;
tunedConfig.det.minSep_mm  = bestMinSep;
tunedConfig.det.roiSize_px = tuner.det.roiSize_px;
tunedConfig.track.maxDisp_mm     = bestMaxDisp;
tunedConfig.track.maxGapFrames   = tuner.track.maxGapFrames;
tunedConfig.track.minTrackLength = bestMinLen;
tunedConfig.track.kalman.processNoise = bestProcessNoise;
tunedConfig.track.kalman.measNoise    = bestMeasNoise;
tunedConfig.motionCorrection     = tuner.motionCorrection;

save(fullfile(tuner.outputFolder, 'tuned_config.mat'), 'tunedConfig');
fprintf('\nSaved to %s/tuned_config.mat\n', tuner.outputFolder);

% Save all sweep data for post-hoc analysis
sweepData.svdStats    = svdStats;
sweepData.detResults  = detResults;
sweepData.trackCounts = trackCounts;
sweepData.trackMedLen = trackMedLen;
sweepData.trackMedSpd = trackMedSpd;
sweepData.tuner       = tuner;
sweepData.bestCutoff  = bestCutoff;
sweepData.bestThreshold = bestThreshold;
sweepData.bestMinSep  = bestMinSep;
sweepData.bestMaxDisp = bestMaxDisp;
sweepData.bestMinLen  = bestMinLen;
sweepData.kalmanCounts  = kalmanCounts;
sweepData.kalmanMedLen  = kalmanMedLen;
sweepData.kalmanMedSpd  = kalmanMedSpd;
sweepData.bestProcessNoise = bestProcessNoise;
sweepData.bestMeasNoise    = bestMeasNoise;
save(fullfile(tuner.outputFolder, 'tuner_sweep_data.mat'), 'sweepData');

fprintf('\n=== Parameter Tuner complete ===\n');
if tuner.useGPU, reset(gpuDevice); end

%% ========================================================================
%  LOCAL HELPERS
%  ========================================================================
function close_if_valid(h)
    if isvalid(h), close(h); end
end
