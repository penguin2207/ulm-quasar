%% QUASAR_PIPELINE.m
% Quantitative Ultrasound Assessment via Sparse Amplitude Recovery
%
% SUSHI-based super-resolution density mapping from VADA plane wave data.
% Processes 1-3 VADA blocks through:
%   1. Beamform to complex IQ (reuses LAT-ULM beamformer)
%   2. SVD clutter filter
%   3. Ensemble power estimation
%   4. Sparse recovery via FISTA (standard SUSHI)
%   5. QUASAR debiased LS amplitude refit (novel contribution)
%   6. Density/concentration map accumulation
%
% The comparison mode runs both L1-only (SUSHI) and L1+LS (QUASAR) on
% every ensemble and reports amplitude statistics to quantify the bias.
%
% Author: Eli Wirth-Apley / Sun Lab, Northeastern University, March 2026

clearvars; close all; clc;

%% ========================================================================
%  USER CONFIGURATION
%  ========================================================================

% --- Input ---
config.mode = 'folder';  % 'single' or 'folder'
config.dataFolder      = 'C:\path\to\VADA_data';                 % <-- set to your .vada data folder
config.baseFilename    = '';
config.modeName        = '.vada';
config.vadaScriptsPath = 'C:\path\to\VevoF2\RFExampleScriptsF2\Matlab';   % <-- VADA SDK scripts (VsiVadaDataRead etc.)
config.outputFolder    = fullfile(config.dataFolder, 'QUASAR_Results');

% --- GPU ---
config.useGPU = true;
config.useCUDA = false;  % true = MEX-CUDA beamformer (requires compiled beamform_pw_das)
config.useRSVD = false;  % true = randomized SVD clutter filter (only when cutoffHigh is empty)

% --- Acquisition ---
% numAngles, numPolarities, and eventsPerFrame are AUTO-DETECTED from data.
% Set config.eventsPerFrame = [] to auto-detect, or override manually.
config.eventsPerFrame = [];         % [] = auto-detect from data

% --- Angle mode ---
config.zeroOnly = true;

% --- Overrides ---
config.sosOverride   = [];
config.pitchOverride = 0.300;

% --- Background ---
config.bgFile   = '';
config.bgFolder = '';

% --- Beamforming ---
% Set to [] for auto-scaling from probe metadata and depth offset.
config.bf.xRange = [];            % [] = auto from RX aperture
config.bf.zRange = [];            % [] = auto from depth offset + data range
config.bf.dx     = [];            % [] = auto ~lambda/5
config.bf.dz     = [];            % [] = auto ~lambda/10

% --- SVD ---
config.svd.cutoffLow  = 2;
config.svd.cutoffHigh = [];

% --- Ensemble ---
config.ensemble.size     = 150;  % Frames per ensemble (50-200 typical)
config.ensemble.overlap  = 0;    % Overlap between ensembles (0 = no overlap)

% --- SUSHI/QUASAR sparse recovery ---
config.sushi.srFactor    = 5;     % Super-resolution factor (pixel subdivision)
config.sushi.lambda      = 0.1;   % L1 regularization (tune per concentration)
config.sushi.maxIter     = 100;   % FISTA iterations
config.sushi.method      = 'fista';
config.sushi.nonNeg      = true;

% --- QUASAR debiased refit ---
config.quasar.enable        = true;
config.quasar.supportThresh = 0;
config.quasar.maxIterCG     = 50;

% --- Decorrelation velocity ---
config.velocity.enable = true;    % Estimate velocity from speckle decorrelation

% --- Comparison mode ---
config.comparison = true;

% --- Cardiac motion correction (optional) ---
% Corrects global rigid tissue motion (cardiac/respiratory) before SVD.
% Uses phase correlation on tissue envelope for sub-pixel shift estimation,
% then Fourier shift theorem to realign complex IQ frames.
config.motionCorrection.enable    = false;       % Set true for in vivo / cardiac
config.motionCorrection.method    = 'phase_corr'; % Phase correlation
config.motionCorrection.refType   = 'rolling';    % 'rolling' or 'first'
config.motionCorrection.refWindow = 10;           % Rolling ref half-width [frames]
config.motionCorrection.maxShift  = 5;            % Max expected shift [pixels]

% --- Per-block diagnostics ---
config.perBlockDiagnostics = true;

% --- Chunk size for beamforming ---
config.chunkSize = 400;

% --- Blanking ---
config.blankSteering = true;
config.blankMargin   = 1.5;
config.blankVoltage  = [];

%% ========================================================================
%  SETUP (reuse LAT-ULM infrastructure)
%  ========================================================================
addpath(genpath(config.vadaScriptsPath));
pipelineDir = fileparts(mfilename('fullpath'));
if isempty(pipelineDir), pipelineDir = pwd; end
addpath(pipelineDir);

if ~exist(config.outputFolder, 'dir'), mkdir(config.outputFolder); end

if config.useGPU
    if gpuDeviceCount == 0
        warning('No GPU found.'); config.useGPU = false;
    else
        g = gpuDevice; reset(g);
        fprintf('GPU: %s (%.1f GB VRAM)\n', g.Name, g.TotalMemory/1e9);
    end
end

fprintf('\n=== QUASAR Pipeline (GPU=%s) ===\n', string(config.useGPU));
pipelineTimer = tic;

%% ========================================================================
%  FILE DISCOVERY AND SERIES SELECTION
%  ========================================================================
switch lower(config.mode)
    case 'single'
        if isempty(config.baseFilename)
            error('config.mode is ''single'' but config.baseFilename is empty.');
        end
        fileList    = {config.baseFilename};
        fileFolders = {config.dataFolder};
        seriesNames = {'single_file'};
        studyNames  = {'single_file'};
        
    case 'folder'
        fprintf('Scanning for .vada files in:\n  %s\n', config.dataFolder);
        vadaFiles = dir(fullfile(config.dataFolder, '**', ['*' config.modeName]));
        if isempty(vadaFiles), error('No %s files found', config.modeName); end
        
        nFound = numel(vadaFiles);
        fileList    = cell(nFound,1);
        fileFolders = cell(nFound,1);
        studyNames  = cell(nFound,1);
        seriesNames = cell(nFound,1);
        fileSizesGB = zeros(nFound,1);
        
        fprintf('Reading metadata from %d file(s)...\n', nFound);
        for i = 1:nFound
            fullName = vadaFiles(i).name;
            fileList{i}    = fullName(1:end-numel(config.modeName));
            fileFolders{i} = vadaFiles(i).folder;
            fileSizesGB(i) = vadaFiles(i).bytes / 1e9;
            
            xmlPath = fullfile(fileFolders{i}, [fileList{i} config.modeName '.xml']);
            params = read_vada_xml_params(xmlPath);
            studyNames{i}  = get_param(params, 'Study-Name', 'UnknownStudy');
            seriesNames{i} = get_param(params, 'Series-Name', 'UnknownSeries');
        end
        
        [~, uIdx] = unique(fileList, 'stable');
        fileList    = fileList(uIdx);
        fileFolders = fileFolders(uIdx);
        studyNames  = studyNames(uIdx);
        seriesNames = seriesNames(uIdx);
        fileSizesGB = fileSizesGB(uIdx);
        nFound = numel(fileList);
        
    otherwise
        error('config.mode must be ''single'' or ''folder''. Got: %s', config.mode);
end

%% Series summary and selection (folder mode)
if strcmpi(config.mode, 'folder')
    seriesKeys = strcat(studyNames, ' | ', seriesNames);
    uniqueSeries = unique(seriesKeys, 'stable');
    nSeries = numel(uniqueSeries);
    estMinPerBlock = 5;  % QUASAR is slower than LAT-ULM
    
    fprintf('\n');
    fprintf('==========================================================================\n');
    fprintf('  SERIES SUMMARY\n');
    fprintf('==========================================================================\n');
    fprintf('  %-4s  %-50s  %6s  %8s  %8s\n', '#', 'Study | Series', 'Blocks', 'Size(GB)', 'Est Time');
    fprintf('  %s\n', repmat('-', 1, 82));
    
    for s = 1:nSeries
        members = find(strcmp(seriesKeys, uniqueSeries{s}));
        nBlk = numel(members);
        totalGB = sum(fileSizesGB(members));
        estMin = nBlk * estMinPerBlock;
        fprintf('  %-4d  %-50s  %6d  %8.1f  %5.0f min\n', ...
            s, uniqueSeries{s}, nBlk, totalGB, estMin);
    end
    
    totalEstMin = nFound * estMinPerBlock;
    fprintf('  %s\n', repmat('-', 1, 82));
    fprintf('  %-4s  %-50s  %6d  %8.1f  %5.0f min\n', ...
        '', 'TOTAL', nFound, sum(fileSizesGB), totalEstMin);
    fprintf('==========================================================================\n\n');
    
    fprintf('Options:\n');
    fprintf('  all     - Process all series\n');
    fprintf('  1,3,5   - Process specific series by number\n');
    fprintf('  1-3     - Process a range of series\n');
    fprintf('  q       - Quit\n\n');
    
    userInput = input('Select series to process: ', 's');
    userInput = strtrim(userInput);
    
    if strcmpi(userInput, 'q'), fprintf('Cancelled.\n'); return; end
    
    if strcmpi(userInput, 'all')
        selectedSeries = 1:nSeries;
    else
        selectedSeries = [];
        parts = strsplit(userInput, ',');
        for p = 1:numel(parts)
            tok = strtrim(parts{p});
            rangeParts = strsplit(tok, '-');
            if numel(rangeParts) == 2
                r1 = str2double(rangeParts{1}); r2 = str2double(rangeParts{2});
                if ~isnan(r1) && ~isnan(r2), selectedSeries = [selectedSeries, r1:r2]; end
            else
                val = str2double(tok);
                if ~isnan(val), selectedSeries = [selectedSeries, val]; end
            end
        end
        selectedSeries = unique(selectedSeries);
        selectedSeries = selectedSeries(selectedSeries >= 1 & selectedSeries <= nSeries);
    end
    
    if isempty(selectedSeries), fprintf('No valid series. Exiting.\n'); return; end
    
    selectedMask = false(nFound, 1);
    for s = selectedSeries
        selectedMask = selectedMask | strcmp(seriesKeys, uniqueSeries{s});
    end
    
    fileList    = fileList(selectedMask);
    fileFolders = fileFolders(selectedMask);
    fileSizesGB = fileSizesGB(selectedMask);
    
    fprintf('\nSelected %d block(s) from %d series:\n', numel(fileList), numel(selectedSeries));
    for s = selectedSeries
        fprintf('  -> %s\n', uniqueSeries{s});
    end
    
    % --- Block-level selection ---
    fprintf('\n  %-4s  %-50s  %8s\n', '#', 'Filename', 'Size(GB)');
    fprintf('  %s\n', repmat('-', 1, 66));
    for i = 1:numel(fileList)
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
                if ~isnan(r1) && ~isnan(r2), selectedBlocks = [selectedBlocks, r1:r2]; end
            else
                val = str2double(tok);
                if ~isnan(val), selectedBlocks = [selectedBlocks, val]; end
            end
        end
        selectedBlocks = unique(selectedBlocks);
        selectedBlocks = selectedBlocks(selectedBlocks >= 1 & selectedBlocks <= numel(fileList));
        
        if isempty(selectedBlocks)
            fprintf('No valid blocks. Exiting.\n'); return;
        end
        
        fileList    = fileList(selectedBlocks);
        fileFolders = fileFolders(selectedBlocks);
        fprintf('\nFiltered to %d block(s).\n', numel(fileList));
    end
end

numBlocks = numel(fileList);
fprintf('\nWill process %d block(s).\n', numBlocks);

%% ========================================================================
%  METADATA + DELAY TABLES
%  ========================================================================
fprintf('\n[Step 0] Loading metadata from first block...\n');

% --- Auto-detect event structure (same logic as LAT_ULM) ---
if isempty(config.eventsPerFrame)
    numProbe = 30;
    [VadaProbe, Param, TxrParam, Config] = VsiVadaDataRead(...
        fileFolders{1}, fileList{1}, 1:numProbe, config.modeName);

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
                    eventsPerFrame = patLen;
                    break;
                end
            end
        end
    end
    if eventsPerFrame == 0
        eventsPerFrame = 6;
        fprintf('  WARNING: Could not auto-detect events/frame. Using default %d.\n', eventsPerFrame);
    end
    config.eventsPerFrame = eventsPerFrame;

    firstFrameAngles = probeAngles(1:eventsPerFrame);
    firstFramePolar  = probePolar(1:eventsPerFrame);
    uAngles = unique(firstFrameAngles, 'stable');
    hasPI = true;
    for a = 1:numel(uAngles)
        if sum(firstFrameAngles == uAngles(a)) < 2
            hasPI = false;
            break;
        end
    end
    config.hasPI = hasPI;
    config.numAngles = numel(uAngles);

    fprintf('  Auto-detected: %d events/frame, %d angles, PI=%s\n', ...
        eventsPerFrame, config.numAngles, string(hasPI));

    VadaTest = VadaProbe(1:eventsPerFrame);
    clear VadaProbe;
else
    [VadaTest, Param, TxrParam, Config] = VsiVadaDataRead(...
        fileFolders{1}, fileList{1}, 1:config.eventsPerFrame, config.modeName);

    probeAngles = zeros(config.eventsPerFrame, 1);
    probePolar  = zeros(config.eventsPerFrame, 1);
    for ev = 1:config.eventsPerFrame
        if isfield(VadaTest(ev).TxDelay, 'angle')
            probeAngles(ev) = VadaTest(ev).TxDelay.angle;
        end
        if isfield(VadaTest(ev).Waveform, 'Channel') && ...
                isfield(VadaTest(ev).Waveform.Channel(1), 'invert')
            probePolar(ev) = VadaTest(ev).Waveform.Channel(1).invert;
        end
    end
    uAngles = unique(probeAngles, 'stable');
    config.numAngles = numel(uAngles);
    hasPI = true;
    for a = 1:numel(uAngles)
        if sum(probeAngles == uAngles(a)) < 2
            hasPI = false; break;
        end
    end
    config.hasPI = hasPI;
    fprintf('  Manual eventsPerFrame=%d, %d angles, PI=%s\n', ...
        config.eventsPerFrame, config.numAngles, string(hasPI));
end

% Pitch
rawPitch = TxrParam.ArrayPitch;
if rawPitch == 0 || isnan(rawPitch), pitch_mm = config.pitchOverride;
elseif rawPitch < 10, pitch_mm = rawPitch;
else, pitch_mm = rawPitch / 1000; end

fs_MHz = Param.SampleFreq;
depthOffset_mm = Param.DepthOffset;
c_xml = Param.SoSMedia;
if ~isempty(config.sosOverride), c = config.sosOverride;
else, c = c_xml; if c==0, c=1540; end; end

fprintf('  %s | pitch=%.4f mm | Fs=%.1f MHz | SoS=%.0f m/s\n', ...
    TxrParam.Name, pitch_mm, fs_MHz, c);

% Detect event structure
angleOrder = zeros(config.eventsPerFrame,1);
polarityOrder = zeros(config.eventsPerFrame,1);
for ev = 1:config.eventsPerFrame
    if isfield(VadaTest(ev).TxDelay,'angle'), angleOrder(ev)=VadaTest(ev).TxDelay.angle; end
    if isfield(VadaTest(ev).Waveform,'Channel') && isfield(VadaTest(ev).Waveform.Channel(1),'invert')
        polarityOrder(ev)=VadaTest(ev).Waveform.Channel(1).invert; end
end
for ev = 1:config.eventsPerFrame
    fprintf('    Evt %d: angle=%+5.1f, inv=%d, Rx=[%d..%d](%d)\n', ev, ...
        angleOrder(ev), polarityOrder(ev), ...
        min(VadaTest(ev).Elements), max(VadaTest(ev).Elements), numel(VadaTest(ev).Elements));
end

uniqueAngles = unique(angleOrder,'stable');
anglePairs = struct('angle',{},'posIdx',{},'negIdx',{},'rxElements',{});
for a = 1:numel(uniqueAngles)
    anglePairs(a).angle = uniqueAngles(a);
    idxs = find(angleOrder==uniqueAngles(a));
    if numel(idxs) >= 2
        if polarityOrder(idxs(1))==0
            anglePairs(a).posIdx=idxs(1); anglePairs(a).negIdx=idxs(2);
        else
            anglePairs(a).posIdx=idxs(2); anglePairs(a).negIdx=idxs(1);
        end
    else
        anglePairs(a).posIdx=idxs(1); anglePairs(a).negIdx=[];
    end
    anglePairs(a).rxElements = VadaTest(anglePairs(a).posIdx).Elements;
end
if hasPI
    config.numPolarities = 2;
else
    config.numPolarities = 1;
end
fprintf('  PI detected: %s, %d angles, %d events/frame\n', string(hasPI), config.numAngles, config.eventsPerFrame);

% --- XML parameter validation ---
% Read all XML parameters and warn if config differs from data
xmlPath = fullfile(fileFolders{1}, [fileList{1} config.modeName '.xml']);
xmlParams = read_vada_xml_params(xmlPath);

fprintf('\n  --- Parameter Validation (config vs XML) ---\n');
warnCount = 0;

% SoS
xmlSoS = str2double(get_param(xmlParams, 'Vada-Mode/Speed-Of-Sound-Media', '0'));
if xmlSoS > 0 && ~isempty(config.sosOverride) && abs(config.sosOverride - xmlSoS) > 1
    fprintf('  WARNING: SoS override (%.0f) differs from XML (%.0f)\n', config.sosOverride, xmlSoS);
    warnCount = warnCount + 1;
end

% Pitch
xmlPitch = get_param(xmlParams, 'Element-Pitch', '');
if ~isempty(xmlPitch)
    xmlPitch_mm = str2double(xmlPitch);
    if abs(xmlPitch_mm - pitch_mm) > 0.001
        fprintf('  WARNING: Pitch used (%.4f mm) differs from XML Element-Pitch (%.4f mm)\n', ...
            pitch_mm, xmlPitch_mm);
        warnCount = warnCount + 1;
    end
end

% Voltage
voltHi = get_param(xmlParams, 'Vada-Mode/Voltage-Rail-High', '');
voltLo = get_param(xmlParams, 'Vada-Mode/Voltage-Rail-Low', '');
txVoltage = [];
if ~isempty(voltHi)
    txVoltage = str2double(voltHi);
    fprintf('  TX voltage: %s-%s%%\n', voltLo, voltHi);
    if ~isempty(config.blankVoltage) && abs(config.blankVoltage - txVoltage) > 1
        fprintf('  WARNING: blankVoltage override (%.0f) differs from XML (%.0f)\n', ...
            config.blankVoltage, txVoltage);
        warnCount = warnCount + 1;
    end
end
config.txVoltage = txVoltage;

% Sequence name
seqName = get_param(xmlParams, 'Vada-Mode/User-Pulse-Sequence-Name', '');
if ~isempty(seqName), fprintf('  Sequence: %s\n', seqName); end

% Number of angles/events
fprintf('  Detected: %d angles, %d events/frame, PI=%s\n', ...
    config.numAngles, config.eventsPerFrame, string(config.hasPI));

% Depth range vs data depth
maxDataDepth_mm = depthOffset_mm + (size(VadaTest(1).Data,1) / (fs_MHz * 2)) * (c * 1e-3);
if ~isempty(config.bf.zRange) && config.bf.zRange(2) > maxDataDepth_mm
    fprintf('  WARNING: config.bf.zRange(2) (%.1f mm) exceeds max data depth (%.1f mm)\n', ...
        config.bf.zRange(2), maxDataDepth_mm);
    warnCount = warnCount + 1;
elseif isempty(config.bf.zRange)
    fprintf('  Depth range: auto (data spans %.1f to %.1f mm)\n', depthOffset_mm, maxDataDepth_mm);
end

if warnCount == 0
    fprintf('  All parameters consistent.\n');
else
    fprintf('  %d warning(s). Review settings above.\n', warnCount);
end
fprintf('  ---\n');

% Zero-only setup
if config.zeroOnly
    zeroAngleIdx = find([anglePairs.angle]==0);
    if isempty(zeroAngleIdx), [~,zeroAngleIdx]=min(abs([anglePairs.angle])); end
    config.zeroAngleIdx = zeroAngleIdx;
    config.zeroEventIdx = anglePairs(zeroAngleIdx).posIdx;
end

% Read TX frequency and data dimensions before clearing VadaTest
txFreq_MHz = 6;
if isfield(VadaTest(1).Waveform, 'Channel') && ~isempty(VadaTest(1).Waveform.Channel)
    if isfield(VadaTest(1).Waveform.Channel(1), 'frequency')
        txFreq_MHz = VadaTest(1).Waveform.Channel(1).frequency;
    end
end
lambda_mm = c * 1e-3 / txFreq_MHz;
nSamplesPerEvent = size(VadaTest(1).Data, 1);
nRxElements = numel(anglePairs(1).rxElements);

clear VadaTest;

fprintf('  TX frequency: %.1f MHz, lambda: %.3f mm\n', txFreq_MHz, lambda_mm);

% Beamforming grid with auto-scaling
elemPos_mm = ((1:TxrParam.ArrayNumElements)-(TxrParam.ArrayNumElements+1)/2)*pitch_mm;

rxElemPos = elemPos_mm(anglePairs(1).rxElements);
rxSpan = max(rxElemPos) - min(rxElemPos);
maxDataDepth_mm = depthOffset_mm + (nSamplesPerEvent / (fs_MHz * 2)) * (c * 1e-3);

if isempty(config.bf.dx),     config.bf.dx = round(lambda_mm / 5, 4);          end
if isempty(config.bf.dz),     config.bf.dz = round(lambda_mm / 10, 4);         end
if isempty(config.bf.xRange)
    config.bf.xRange = [min(rxElemPos) - rxSpan*0.2, max(rxElemPos) + rxSpan*0.2];
end
if isempty(config.bf.zRange)
    config.bf.zRange = [max(depthOffset_mm, 0.5), maxDataDepth_mm];
end

fprintf('  Grid: dx=%.4f mm, dz=%.4f mm\n', config.bf.dx, config.bf.dz);
fprintf('  X: [%.1f, %.1f] mm, Z: [%.1f, %.1f] mm\n', ...
    config.bf.xRange(1), config.bf.xRange(2), config.bf.zRange(1), config.bf.zRange(2));

xGrid = config.bf.xRange(1):config.bf.dx:config.bf.xRange(2);
zGrid = config.bf.zRange(1):config.bf.dz:config.bf.zRange(2);
nX = numel(xGrid); nZ = numel(zGrid);

if config.zeroOnly
    rxPos_mm = elemPos_mm(anglePairs(config.zeroAngleIdx).rxElements);
    delayTables = {beamform_planewave_gpu([], rxPos_mm, anglePairs(config.zeroAngleIdx).angle, ...
        xGrid, zGrid, fs_MHz, c, depthOffset_mm, [])};
    nRx = numel(anglePairs(config.zeroAngleIdx).rxElements);
else
    delayTables = cell(config.numAngles,1);
    for a = 1:config.numAngles
        rxPos_mm = elemPos_mm(anglePairs(a).rxElements);
        delayTables{a} = beamform_planewave_gpu([], rxPos_mm, anglePairs(a).angle, ...
            xGrid, zGrid, fs_MHz, c, depthOffset_mm, []);
    end
    nRx = numel(anglePairs(1).rxElements);
end

fprintf('  BF grid: %d x %d = %d pixels\n', nX, nZ, nX*nZ);

%% ========================================================================
%  BACKGROUND (same as LAT-ULM)
%  ========================================================================
bgMeanIQ = [];
if ~isempty(config.bgFile)
    fprintf('\n[Step 1b] Computing background...\n');
    bgFolder = config.bgFolder; if isempty(bgFolder), bgFolder = config.dataFolder; end
    try
        [VadaBg,~,~,BgConfig] = VsiVadaDataRead(bgFolder, config.bgFile, ...
            1:config.eventsPerFrame, config.modeName);
        numBgEvents = numel(BgConfig.PulseSequences(1).Events);
        numBgFrames = min(floor(numBgEvents/config.eventsPerFrame), 500);
        clear VadaBg;
        [VadaBg,~,~,~] = VsiVadaDataRead(bgFolder, config.bgFile, ...
            1:numBgFrames*config.eventsPerFrame, config.modeName);

        % Select beamformer for background stack
        if isfield(config, 'useCUDA') && config.useCUDA
            bg_bf_fn = @beamform_cuda;
        else
            bg_bf_fn = @beamform_planewave_gpu;
        end

        % Beamform into full IQ stack (needed for motion correction)
        bgStack = zeros(nZ, nX, numBgFrames, 'single');
        for iFrame = 1:numBgFrames
            baseEvt = (iFrame-1)*config.eventsPerFrame;
            if config.zeroOnly
                evIdx = baseEvt + config.zeroEventIdx;
                rfData = single(VadaBg(evIdx).Data);
                rxPos = elemPos_mm(anglePairs(config.zeroAngleIdx).rxElements);
                bfImg = bg_bf_fn(rfData, rxPos, ...
                    anglePairs(config.zeroAngleIdx).angle, xGrid, zGrid, ...
                    fs_MHz, c, depthOffset_mm, delayTables{1});
            else
                bfImg = complex(zeros(nZ, nX, 'single'));
                for a = 1:config.numAngles
                    rfPos = single(VadaBg(baseEvt + anglePairs(a).posIdx).Data);
                    if config.hasPI
                        rfNeg = single(VadaBg(baseEvt + anglePairs(a).negIdx).Data);
                        rfData = rfPos + rfNeg;
                    else
                        rfData = rfPos;
                    end
                    rxPos = elemPos_mm(anglePairs(a).rxElements);
                    bfA = bg_bf_fn(rfData, rxPos, anglePairs(a).angle, ...
                        xGrid, zGrid, fs_MHz, c, depthOffset_mm, delayTables{a});
                    bfImg = bfImg + single(bfA);
                end
            end
            bgStack(:,:,iFrame) = single(bfImg);
        end
        clear VadaBg;

        % Motion-correct background before averaging (in vivo: tissue moves)
        if isfield(config, 'motionCorrection') && config.motionCorrection.enable
            fprintf('  Motion-correcting background...\n');
            mcP.method    = config.motionCorrection.method;
            mcP.refType   = config.motionCorrection.refType;
            mcP.refWindow = config.motionCorrection.refWindow;
            mcP.maxShift  = config.motionCorrection.maxShift;
            [bgShifts, bgMcDiag] = estimate_tissue_motion(bgStack, mcP);
            bgStack = apply_motion_correction(bgStack, bgShifts, config.useGPU);
            fprintf('  BG motion: max=%.2f px, mean=%.2f px\n', ...
                bgMcDiag.maxDisp_px, bgMcDiag.meanDisp_px);
        end

        bgMeanIQ = mean(bgStack, 3);
        clear bgStack;
        fprintf('  Background: %d frames (motion-corrected=%s)\n', numBgFrames, ...
            string(isfield(config, 'motionCorrection') && config.motionCorrection.enable));
    catch ME
        fprintf('  Background failed: %s\n', ME.message);
    end
end

%% ========================================================================
%  BUILD PSF MODEL
%  ========================================================================
fprintf('\n[Step 2] Building PSF model...\n');

srFactor = config.sushi.srFactor;
srPixel_mm = config.bf.dx / srFactor;  % Use lateral pixel as reference
srNZ = nZ * srFactor;
srNX = nX * srFactor;

[psf, psfParams] = build_sushi_psf(txFreq_MHz, c, pitch_mm, nRx, ...
    srPixel_mm, [srNZ, srNX], 0);

fprintf('  PSF: axial FWHM=%.3f mm (%.1f SR px), lateral FWHM=%.3f mm (%.1f SR px)\n', ...
    psfParams.axialFWHM_mm, psfParams.axialFWHM_px, ...
    psfParams.lateralFWHM_mm, psfParams.lateralFWHM_px);
fprintf('  SR grid: %d x %d (%.3f mm pixel, %dx factor)\n', srNZ, srNX, srPixel_mm, srFactor);

% SR coordinate vectors (needed for per-block diagnostics and final output)
srX = linspace(config.bf.xRange(1), config.bf.xRange(2), srNX);
srZ = linspace(config.bf.zRange(1), config.bf.zRange(2), srNZ);

%% ========================================================================
%  PROCESS BLOCKS -> ENSEMBLES -> SPARSE RECOVERY
%  ========================================================================
fprintf('\n[Step 3] Processing blocks...\n');

% Compute effective frame rate from event timing
try
    [VadaTmp,~,~,~] = VsiVadaDataRead(fileFolders{1}, fileList{1}, ...
        1:(2*config.eventsPerFrame), config.modeName);
    framePeriod_ms = VadaTmp(config.eventsPerFrame+1).Timestamp - VadaTmp(1).Timestamp;
    frameRate = 1000 / framePeriod_ms;
    clear VadaTmp;
catch
    frameRate = 1000 / (config.eventsPerFrame * 0.150);  % Fallback estimate
end
fprintf('  Effective frame rate: %.0f Hz\n', frameRate);

ensembleSize = config.ensemble.size;
ensembleOverlap = config.ensemble.overlap;
ensembleStep = ensembleSize - ensembleOverlap;

% Accumulation maps
sushiDensity   = zeros(srNZ, srNX, 'single');  % L1-only (SUSHI)
quasarDensity  = zeros(srNZ, srNX, 'single');  % L1+LS (QUASAR)
velAccum       = zeros(srNZ, srNX, 'single');  % Velocity accumulator
velCount       = zeros(srNZ, srNX, 'single');  % Valid velocity pixel count
nEnsembles     = 0;

% Comparison stats
compStats_cells = {};

fistaOpts.maxIter = config.sushi.maxIter;
fistaOpts.nonNeg  = config.sushi.nonNeg;
fistaOpts.useGPU  = config.useGPU;
fistaOpts.verbose = false;

quasarOpts.maxIterCG     = config.quasar.maxIterCG;
quasarOpts.supportThresh = config.quasar.supportThresh;
quasarOpts.useGPU        = config.useGPU;

for iBlock = 1:numBlocks
    fprintf('\n  --- Block %d/%d: %s ---\n', iBlock, numBlocks, fileList{iBlock});
    
    % Per-block accumulators for diagnostics
    blockSushi  = zeros(srNZ, srNX, 'single');
    blockQuasar = zeros(srNZ, srNX, 'single');
    blockVelA   = zeros(srNZ, srNX, 'single');
    blockVelC   = zeros(srNZ, srNX, 'single');
    nBlockEns   = 0;
    
    % Load metadata for this block
    [~,~,~,BlockConfig] = VsiVadaDataRead(fileFolders{iBlock}, fileList{iBlock}, ...
        1:config.eventsPerFrame, config.modeName);
    numTotalEvents = numel(BlockConfig.PulseSequences(1).Events);
    numCompoundFrames = floor(numTotalEvents / config.eventsPerFrame);
    numBlockEnsembles = floor((numCompoundFrames - ensembleSize) / ensembleStep) + 1;
    
    fprintf('    %d frames, %d ensembles (size=%d, step=%d)\n', ...
        numCompoundFrames, numBlockEnsembles, ensembleSize, ensembleStep);
    
    % Process in chunks (beamform), then extract ensembles
    % Load and beamform entire block chunk-by-chunk
    numChunks = ceil(numCompoundFrames / config.chunkSize);
    IQ_full = zeros(nZ, nX, numCompoundFrames, 'single');

    % Select beamformer: CUDA MEX kernel or MATLAB gpuArray
    if isfield(config, 'useCUDA') && config.useCUDA
        bf_fn = @beamform_cuda;
    else
        bf_fn = @beamform_planewave_gpu;
    end

    for iChunk = 1:numChunks
        frameStart = (iChunk-1)*config.chunkSize + 1;
        frameEnd = min(iChunk*config.chunkSize, numCompoundFrames);
        nFrames = frameEnd - frameStart + 1;
        eventStart = (frameStart-1)*config.eventsPerFrame + 1;
        eventEnd = frameEnd * config.eventsPerFrame;
        
        try
            [VadaChunk,~,~,~] = VsiVadaDataRead(fileFolders{iBlock}, fileList{iBlock}, ...
                eventStart:eventEnd, config.modeName);
        catch ME
            fprintf('      Chunk %d skipped: %s\n', iChunk, ME.message);
            continue;
        end
        
        for iFrame = 1:nFrames
            baseEvt = (iFrame-1)*config.eventsPerFrame;
            if config.zeroOnly
                evIdx = baseEvt + config.zeroEventIdx;
                rfData = single(VadaChunk(evIdx).Data);
                rxPos = elemPos_mm(anglePairs(config.zeroAngleIdx).rxElements);
                bfImg = bf_fn(rfData, rxPos, ...
                    anglePairs(config.zeroAngleIdx).angle, ...
                    xGrid, zGrid, fs_MHz, c, depthOffset_mm, delayTables{1});
                IQ_full(:,:,frameStart+iFrame-1) = single(bfImg);
            else
                % All angles with coherent compound (PI if available)
                compImg = complex(zeros(nZ, nX, 'single'));
                for a = 1:config.numAngles
                    rfPos = single(VadaChunk(baseEvt + anglePairs(a).posIdx).Data);
                    if config.hasPI
                        rfNeg = single(VadaChunk(baseEvt + anglePairs(a).negIdx).Data);
                        if config.blankSteering && anglePairs(a).angle ~= 0
                            blankInfo = compute_steering_blanking(anglePairs(a).angle, ...
                                numel(anglePairs(a).rxElements), ...
                                elemPos_mm(2)-elemPos_mm(1), c, fs_MHz, ...
                                config.blankMargin, config.txVoltage);
                            if blankInfo.nBlank > 0 && blankInfo.nBlank < size(rfPos,1)
                                rfPos(1:blankInfo.nBlank, :) = 0;
                                rfNeg(1:blankInfo.nBlank, :) = 0;
                            end
                        end
                        rfData = rfPos + rfNeg;
                    else
                        rfData = rfPos;
                    end
                    rxPos = elemPos_mm(anglePairs(a).rxElements);
                    bfImg = bf_fn(rfData, rxPos, anglePairs(a).angle, ...
                        xGrid, zGrid, fs_MHz, c, depthOffset_mm, delayTables{a});
                    compImg = compImg + single(bfImg);
                end
                IQ_full(:,:,frameStart+iFrame-1) = compImg;
            end
        end
        clear VadaChunk;
        
        if mod(iChunk, 5) == 0
            fprintf('      Beamformed chunk %d/%d\n', iChunk, numChunks);
        end
    end
    
    % Motion correction (if enabled) — must precede SVD for clean clutter separation
    if isfield(config, 'motionCorrection') && config.motionCorrection.enable
        fprintf('    Estimating tissue motion...\n');
        mcParams.method    = config.motionCorrection.method;
        mcParams.refType   = config.motionCorrection.refType;
        mcParams.refWindow = config.motionCorrection.refWindow;
        mcParams.maxShift  = config.motionCorrection.maxShift;

        [mcShifts, mcDiag] = estimate_tissue_motion(IQ_full, mcParams);
        IQ_full = apply_motion_correction(IQ_full, mcShifts, config.useGPU);

        fprintf('    Motion correction: max=%.2f px, mean=%.2f px\n', ...
            mcDiag.maxDisp_px, mcDiag.meanDisp_px);
    end

    % Background subtraction
    if ~isempty(bgMeanIQ)
        for f = 1:numCompoundFrames
            IQ_full(:,:,f) = IQ_full(:,:,f) - bgMeanIQ;
        end
    end

    % Process ensembles (SVD per-ensemble to avoid GPU memory overflow)
    % Full-block SVD would require ~nZ*nX*nFrames*8 bytes on GPU.
    % Per-ensemble SVD on 150 frames is ~500 MB, well within GPU limits.
    fprintf('    Processing ensembles (SVD per-ensemble)');
    for iEns = 1:numBlockEnsembles
        ensStart = (iEns-1)*ensembleStep + 1;
        ensEnd = ensStart + ensembleSize - 1;
        if ensEnd > numCompoundFrames, break; end
        
        % SVD clutter filter on this ensemble only
        ensData = svd_clutter_filter_gpu(IQ_full(:,:,ensStart:ensEnd), ...
            config.svd.cutoffLow, config.svd.cutoffHigh, config.useGPU);
        
        % Decorrelation velocity (before time-averaging destroys temporal info)
        if config.velocity.enable
            [ensVel, ~] = estimate_decorrelation_velocity(ensData, ...
                frameRate, psfParams.lambda_mm, srFactor, config.useGPU);
            validMask = ensVel > 0 & ~isnan(ensVel);
            velAccum(validMask) = velAccum(validMask) + ensVel(validMask);
            velCount(validMask) = velCount(validMask) + 1;
            blockVelA(validMask) = blockVelA(validMask) + ensVel(validMask);
            blockVelC(validMask) = blockVelC(validMask) + 1;
        end
        
        % Compute ensemble power image (mean of |IQ|^2)
        powerImg = mean(abs(ensData).^2, 3);
        clear ensData;
        
        % Upsample to SR grid (bicubic interpolation)
        powerSR = imresize(powerImg, [srNZ, srNX], 'bicubic');
        powerSR = max(powerSR, 0);  % Enforce non-negativity
        
        % Normalize for stable FISTA convergence
        normFactor = max(powerSR(:)) + eps;
        powerNorm = powerSR / normFactor;
        
        % --- FISTA sparse recovery (standard SUSHI) ---
        [x_fista, fistaInfo] = sushi_sparse_recovery(powerNorm, psf, ...
            config.sushi.lambda, config.sushi.method, fistaOpts);
        
        % Rescale
        x_fista = x_fista * normFactor;
        
        % Accumulate SUSHI density (global + per-block)
        sushiDensity = sushiDensity + x_fista;
        blockSushi = blockSushi + x_fista;
        
        % --- QUASAR debiased refit ---
        if config.quasar.enable
            [x_quasar, quasarInfo] = quasar_refit(powerNorm, psf, ...
                x_fista/normFactor, quasarOpts);
            x_quasar = x_quasar * normFactor;
            quasarDensity = quasarDensity + x_quasar;
            blockQuasar = blockQuasar + x_quasar;
        end
        
        nEnsembles = nEnsembles + 1;
        nBlockEns = nBlockEns + 1;
        
        % Store comparison stats
        if config.comparison
            cs.fistaSupport = fistaInfo.support;
            cs.fistaMaxAmp = fistaInfo.maxVal * normFactor;
            cs.fistaTime = fistaInfo.time;
            cs.block = iBlock;
            cs.ensemble = iEns;
            if config.quasar.enable
                cs.quasarAmpRatio = quasarInfo.amplitudeRatio;
                cs.quasarResidual = quasarInfo.residualNorm;
                cs.quasarTime = quasarInfo.time;
            else
                cs.quasarAmpRatio = NaN;
                cs.quasarResidual = NaN;
                cs.quasarTime = NaN;
            end
            compStats_cells{end+1} = cs; %#ok
        end
        
        if mod(iEns, 10) == 0, fprintf('.'); end
    end
    fprintf(' done (%d ensembles)\n', min(iEns, numBlockEnsembles));
    
    % Per-block diagnostics
    if config.perBlockDiagnostics
        blockVelMap = blockVelA ./ (blockVelC + eps);
        blockVelMap(blockVelC == 0) = 0;
        save_block_diagnostics(config.outputFolder, fileList{iBlock}, iBlock, ...
            [], config, xGrid, zGrid, fs_MHz, c, frameRate, 'quasar', ...
            blockSushi, blockQuasar, blockVelMap, nBlockEns, srX, srZ, []);
    end
    
    clear IQ_full;
end

fprintf('\n  Total ensembles processed: %d\n', nEnsembles);

%% ========================================================================
%  COMPARISON STATISTICS (SUSHI vs QUASAR)
%  ========================================================================
if ~isempty(compStats_cells)
    compStats = [compStats_cells{:}];
else
    compStats = [];
end

if config.comparison && ~isempty(compStats)
    fprintf('\n[Step 4] SUSHI vs QUASAR comparison...\n');
    
    ampRatios = [compStats.quasarAmpRatio];
    fistaSupports = [compStats.fistaSupport];
    fistaTimes = [compStats.fistaTime];
    quasarTimes = [compStats.quasarTime];
    
    fprintf('  %-30s %-15s %-15s\n', 'Metric', 'Mean', 'Std');
    fprintf('  %s\n', repmat('-', 1, 60));
    fprintf('  %-30s %-15.1f %-15.1f\n', 'FISTA support (pixels)', mean(fistaSupports), std(fistaSupports));
    fprintf('  %-30s %-15.2f %-15.2f\n', 'QUASAR amp ratio', nanmean(ampRatios), nanstd(ampRatios));
    fprintf('  %-30s %-15.3f %-15.3f\n', 'FISTA time (s)', mean(fistaTimes), std(fistaTimes));
    fprintf('  %-30s %-15.3f %-15.3f\n', 'QUASAR refit time (s)', nanmean(quasarTimes), nanstd(quasarTimes));
    
    fprintf('\n  Amplitude ratio interpretation:\n');
    fprintf('    ratio > 1.0 : FISTA underestimates amplitudes (L1 bias confirmed)\n');
    fprintf('    ratio ~ 1.0 : minimal bias (lambda may be too low)\n');
    fprintf('    ratio >> 2.0: significant bias correction by QUASAR\n');
    
    % Plot comparison
    figure('Name', 'SUSHI vs QUASAR', 'Position', [100 100 1200 400]);
    
    subplot(1,3,1);
    histogram(ampRatios, 20, 'FaceColor', [0.2 0.6 0.3]);
    xlabel('QUASAR / FISTA Amplitude Ratio');
    ylabel('Ensemble Count');
    title(sprintf('Amplitude Bias (mean=%.2f)', nanmean(ampRatios)));
    xline(1, 'r--', 'No bias', 'LineWidth', 1.5);
    
    subplot(1,3,2);
    histogram(fistaSupports, 20, 'FaceColor', [0.3 0.4 0.8]);
    xlabel('Support Size (pixels)');
    ylabel('Ensemble Count');
    title(sprintf('FISTA Support (mean=%.0f)', mean(fistaSupports)));
    
    subplot(1,3,3);
    scatter(fistaTimes, quasarTimes, 20, 'filled');
    xlabel('FISTA Time [s]'); ylabel('QUASAR Refit Time [s]');
    title('Computation Time'); grid on;
    hold on; plot([0 max(fistaTimes)], [0 max(fistaTimes)], 'k--');
    
    sgtitle('SUSHI (L1) vs QUASAR (L1 + LS Refit)');
    saveas(gcf, fullfile(config.outputFolder, 'sushi_vs_quasar_comparison.png'));
end

%% ========================================================================
%  VISUALIZATION
%  ========================================================================
fprintf('\n[Step 5] Visualization...\n');

% srX and srZ already defined before block loop

% Compute mean velocity map
velMap = velAccum ./ (velCount + eps);
velMap(velCount == 0) = 0;

% --- Figure 1: Density comparison ---
figure('Name', 'QUASAR Density Maps', 'Position', [50 50 1400 500], 'Visible', 'off');

subplot(1,3,1);
imagesc(srX, srZ, log10(sushiDensity / nEnsembles + 1));
axis image; colormap(gca, hot); colorbar;
xlabel('Lateral [mm]'); ylabel('Axial [mm]');
title(sprintf('SUSHI (L1 only)\n%d ensembles, \\lambda=%.3f', nEnsembles, config.sushi.lambda));

if config.quasar.enable
    subplot(1,3,2);
    imagesc(srX, srZ, log10(quasarDensity / nEnsembles + 1));
    axis image; colormap(gca, hot); colorbar;
    xlabel('Lateral [mm]'); ylabel('Axial [mm]');
    title(sprintf('QUASAR (L1 + LS refit)\n%d ensembles', nEnsembles));
    
    subplot(1,3,3);
    ratioMap = quasarDensity ./ (sushiDensity + eps);
    ratioMap(sushiDensity == 0) = 0;
    imagesc(srX, srZ, ratioMap);
    axis image; colorbar; caxis([0 5]);
    xlabel('Lateral [mm]'); ylabel('Axial [mm]');
    title(sprintf('QUASAR/SUSHI Ratio (mean=%.2f)', mean(ratioMap(sushiDensity>0))));
    colormap(gca, 'jet');
end
sgtitle('QUASAR Density');
saveas(gcf, fullfile(config.outputFolder, 'quasar_density_maps.png'));

% --- Figure 2: Velocity map ---
if config.velocity.enable
    figure('Name', 'Velocity Map', 'Position', [50 50 900 500], 'Visible', 'off');
    
    subplot(1,2,1);
    imagesc(srX, srZ, velMap);
    axis image; colorbar; caxis([0 30]);
    colormap(gca, jet);
    xlabel('Lateral [mm]'); ylabel('Axial [mm]');
    title(sprintf('Decorrelation Velocity [mm/s]\n%d ensembles', nEnsembles));
    
    subplot(1,2,2);
    % Overlay velocity on density (mask to signal region)
    if config.quasar.enable
        densNorm = quasarDensity / (max(quasarDensity(:)) + eps);
    else
        densNorm = sushiDensity / (max(sushiDensity(:)) + eps);
    end
    velMasked = velMap;
    velMasked(densNorm < 0.1) = NaN;  % Only show velocity where there's signal
    imagesc(srX, srZ, velMasked);
    axis image; colorbar; caxis([0 30]);
    colormap(gca, jet);
    xlabel('Lateral [mm]'); ylabel('Axial [mm]');
    title('Velocity (density-masked)');
    
    sgtitle('QUASAR Velocity Estimation');
    saveas(gcf, fullfile(config.outputFolder, 'quasar_velocity_map.png'));
end

% --- Figure 3: Statistics summary ---
figure('Name', 'Statistics', 'Position', [50 50 1200 400], 'Visible', 'off');

subplot(1,3,1);
if config.quasar.enable
    densForProfile = quasarDensity / nEnsembles;
else
    densForProfile = sushiDensity / nEnsembles;
end
midZ = round(size(densForProfile,1)/2);
zBand = max(1,midZ-20):min(size(densForProfile,1),midZ+20);
latProfile = mean(densForProfile(zBand,:), 1);
latProfile = latProfile / (max(latProfile) + eps);
plot(srX, latProfile, 'b-', 'LineWidth', 1.5);
xlabel('Lateral [mm]'); ylabel('Normalized Density');
title('Lateral Profile (mid-depth)'); grid on;
yline(0.5, 'r--', 'FWHM');

subplot(1,3,2);
if config.velocity.enable && any(velMap(:) > 0)
    velValid = velMap(velMap > 0 & velMap < 100);
    histogram(velValid(:), 50, 'FaceColor', [0.3 0.5 0.8]);
    xlabel('Velocity [mm/s]'); ylabel('Count');
    title(sprintf('Velocity Distribution\nmedian=%.1f mm/s', median(velValid)));
else
    text(0.5, 0.5, 'No velocity data', 'HorizontalAlignment', 'center');
    axis off;
end

subplot(1,3,3);
if ~isempty(compStats)
    ampRatios_plot = [compStats.quasarAmpRatio];
    ampRatios_plot = ampRatios_plot(~isnan(ampRatios_plot));
    histogram(ampRatios_plot, 20, 'FaceColor', [0.2 0.6 0.3]);
    xlabel('QUASAR/SUSHI Ratio'); ylabel('Count');
    title(sprintf('Amplitude Bias\nmean=%.2f', mean(ampRatios_plot)));
    xline(1, 'r--', 'No bias');
else
    text(0.5, 0.5, 'No comparison data', 'HorizontalAlignment', 'center');
    axis off;
end

sgtitle(sprintf('QUASAR Statistics: %d blocks, %d ensembles', numBlocks, nEnsembles));
saveas(gcf, fullfile(config.outputFolder, 'quasar_statistics.png'));
close all;

%% ========================================================================
%  SAVE
%  ========================================================================
fprintf('\n[Saving]...\n');

R.config          = config;
R.sushiDensity    = sushiDensity;
R.quasarDensity   = quasarDensity;
R.velocityMap     = velMap;
R.velAccum        = velAccum;
R.velCount        = velCount;
R.nEnsembles      = nEnsembles;
R.frameRate       = frameRate;
R.compStats       = compStats;
R.srX             = srX;
R.srZ             = srZ;
R.srPixel_mm      = srPixel_mm;
R.psfParams       = psfParams;
R.Param           = Param;
R.TxrParam        = TxrParam;
R.processingTime  = toc(pipelineTimer);

save(fullfile(config.outputFolder, 'QUASAR_results.mat'), 'R', '-v7.3');

fprintf('\nSaved to %s\n', config.outputFolder);
fprintf('=== QUASAR complete: %d blocks, %d ensembles in %.1f min ===\n', ...
    numBlocks, nEnsembles, toc(pipelineTimer)/60);

%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%  ========================================================================

function params = read_vada_xml_params(xmlPath)
params = struct();
if ~exist(xmlPath, 'file'), return; end
text = fileread(xmlPath);
tokens = regexp(text, '<parameter\s+name="([^"]+)"\s+value="([^"]*)"', 'tokens');
for i = 1:numel(tokens)
    fn = strrep(strrep(tokens{i}{1}, '-', '_'), '/', '_');
    params.(fn) = tokens{i}{2};
end
end

function val = get_param(params, name, default)
fn = strrep(strrep(name, '-', '_'), '/', '_');
if isstruct(params) && isfield(params, fn)
    val = params.(fn);
else
    val = default;
end
end
