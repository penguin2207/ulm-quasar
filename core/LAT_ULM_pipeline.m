%% LAT_ULM_PIPELINE.m
% Complete LAT-ULM pipeline for VEVO F2 / VADA data.
% GPU-accelerated (NVIDIA A6000). Supports single block or batch folder.
%
% MODES:
%   'single' - Process one specified VADA file
%   'folder' - Auto-discover and process ALL .vada blocks in a folder
%
% Tracking runs per-block (68s inter-block gaps prevent cross-block linking).
% Localizations and tracks accumulate across blocks for final SR rendering.
%
% References:
%   Hingot et al., Sci Rep 2019;9:2456
%   Heiles et al. (PALA), Nat Biomed Eng 2022;6(5):605-616
%   Demene et al., IEEE TMI 2015;34(11):2271-2285
%   Tang et al., IEEE TUFFC 2020;67(9):1738-1751
%
% Author: Eli Wirth-Apley / Sun Lab, Northeastern University, February 2026

clearvars; close all; clc;

% Global config override hook: wrapper scripts (e.g., ULM3_LAT_ULM_QC_RUNNER)
% can set this global struct before calling `run('LAT_ULM_pipeline.m')` to
% patch specific config fields without hand-editing this file.
% `clearvars` (line above) doesn't touch globals, but the global must be
% re-declared here to be accessible after clearvars.
global LAT_ULM_CONFIG_OVERRIDE %#ok<GVMIS>

%% ========================================================================
%  USER CONFIGURATION
%  ========================================================================

% --- Input mode ---
% 'single' : process one file specified by config.baseFilename
% 'folder' : auto-discover all .vada files in config.dataFolder
config.mode = 'folder';  % <-- CHANGE: 'single' or 'folder'

% --- Paths ---
config.dataFolder      = 'C:\path\to\VADA_data';                 % <-- set to your .vada data folder
config.baseFilename    = '';  % Only used in 'single' mode (no extension)
config.modeName        = '.vada';
config.vadaScriptsPath = 'C:\path\to\VevoF2\RFExampleScriptsF2\Matlab';   % <-- VADA SDK scripts (VsiVadaDataRead etc.)
config.outputFolder    = fullfile(config.dataFolder, 'LAT_ULM_Results');

% --- GPU ---
config.useGPU = true;
config.useCUDA = false;  % true = MEX-CUDA beamformer (requires compiled beamform_pw_das)
config.useRSVD = false;  % true = randomized SVD clutter filter (only when cutoffHigh is empty)

% --- Acquisition ---
% numAngles and eventsPerFrame are AUTO-DETECTED from the first block.
% Set config.eventsPerFrame = [] to auto-detect, or override manually.
config.eventsPerFrame = [];         % [] = auto-detect from data
config.chunkSize     = 400;        % Compound frames per processing chunk

% --- Angle mode ---
% false = use all angles with PI + coherent compound
% true  = use only 0deg event, no PI, complex IQ for SVD
config.zeroOnly = true;

% --- Steering RF blanking (for multi-angle mode only) ---
% At steered angles, TX cross-talk contaminates the first N RF samples
% asymmetrically for +/- PI pulses. Blanking these samples before PI
% summation improves cancellation quality at steered angles.
% Margin auto-scales with TX voltage (higher power = more nonlinearity).
config.blankSteering = true;     % Enable RF blanking at steered angles
config.blankMargin   = 1.5;      % Base safety margin (auto-scaled by voltage)
config.blankVoltage  = [];       % TX voltage override [%]. [] = read from XML.

% --- Background subtraction ---
% Path to a VADA block acquired with water/PBS only (no bubbles).
% The mean complex IQ frame is subtracted from every bubble frame before SVD.
% Set to '' to disable background subtraction.
config.bgFile   = '';  % e.g., 'water_baseline_block'  (no extension)
config.bgFolder = '';  % e.g., same as dataFolder. Leave '' to use config.dataFolder

% --- Block quality control ---
% Auto-screen blocks for bubble presence (informational, no filtering).
% Results saved in results.blockQC for post-hoc review.

% --- Overrides (applied when XML reports 0 or wrong values) ---
% Speed of sound: XML may report the F2 system setting (e.g., 1480 m/s)
% rather than the actual medium. Set to [] to use XML value, or override:
config.sosOverride   = 1540;       % [m/s] Use 1540 for agarose, [] for XML value
config.pitchOverride = 0.300;      % [mm] Used if ArrayPitch=0 in metadata

% --- Beamforming ---
% Set to [] for auto-scaling from probe metadata and depth offset.
% Auto-scaling uses TX frequency, RX aperture, and depth offset to compute
% appropriate grid ranges and pixel sizes for any probe.
config.bf.xRange = [];            % Lateral [mm]. [] = auto from RX aperture
config.bf.zRange = [];            % Axial [mm].   [] = auto from depth offset + data range
config.bf.dx     = [];            % Lateral pixel [mm]. [] = auto ~lambda/5
config.bf.dz     = [];            % Axial pixel [mm].   [] = auto ~lambda/10

% --- SVD clutter filter ---
config.svd.cutoffLow  = 5;
config.svd.cutoffHigh = [];

% --- Detection ---
config.det.method     = 'threshold';
config.det.threshold  = 5;
config.det.minSep_mm  = [];       % [] = auto ~1 wavelength
config.det.roiSize_px = 7;

% --- Tracking ---
config.track.maxDisp_mm     = 0.500;
config.track.maxGapFrames   = 3;
config.track.minTrackLength = 5;
config.track.kalman.processNoise = 0.01;
config.track.kalman.measNoise    = 0.05;

% --- Super-resolution ---
config.sr.pixelSize_um = [];      % [] = auto ~lambda/5

% --- Cardiac motion correction (optional) ---
% Corrects global rigid tissue motion (cardiac/respiratory) before SVD.
% Uses phase correlation on tissue envelope for sub-pixel shift estimation,
% then Fourier shift theorem to realign complex IQ frames.
% Reference: Demene et al., IEEE TMI 2015;34(11):2271-2285
config.motionCorrection.enable    = false;       % Set true for in vivo / cardiac
config.motionCorrection.method    = 'phase_corr'; % Phase correlation
config.motionCorrection.refType   = 'rolling';    % 'rolling' or 'first'
config.motionCorrection.refWindow = 10;           % Rolling ref half-width [frames]
config.motionCorrection.maxShift  = 5;            % Max expected shift [pixels]

% --- Per-block diagnostics ---
% Save a diagnostic PNG and summary per block to outputDir/blocks/
config.perBlockDiagnostics = true;

%% ========================================================================
%  APPLY WRAPPER CONFIG OVERRIDES (if any)
%  ========================================================================
if ~isempty(LAT_ULM_CONFIG_OVERRIDE) && isstruct(LAT_ULM_CONFIG_OVERRIDE)
    fn = fieldnames(LAT_ULM_CONFIG_OVERRIDE);
    fprintf('\n*** LAT_ULM_CONFIG_OVERRIDE applied (%d fields) ***\n', numel(fn));
    for iOv = 1:numel(fn)
        val = LAT_ULM_CONFIG_OVERRIDE.(fn{iOv});
        if isstruct(val) && isfield(config, fn{iOv}) && isstruct(config.(fn{iOv}))
            % Deep merge: patch fields onto existing substruct (preserves defaults)
            sub_fn = fieldnames(val);
            for jOv = 1:numel(sub_fn)
                sv = val.(sub_fn{jOv});
                if isstruct(sv) && isfield(config.(fn{iOv}), sub_fn{jOv}) && isstruct(config.(fn{iOv}).(sub_fn{jOv}))
                    % One more level deep (e.g., track.kalman)
                    ssf = fieldnames(sv);
                    for kk = 1:numel(ssf)
                        config.(fn{iOv}).(sub_fn{jOv}).(ssf{kk}) = sv.(ssf{kk});
                        fprintf('  config.%s.%s.%s = %g\n', fn{iOv}, sub_fn{jOv}, ssf{kk}, sv.(ssf{kk}));
                    end
                else
                    config.(fn{iOv}).(sub_fn{jOv}) = sv;
                    if isnumeric(sv)
                        fprintf('  config.%s.%s = %s\n', fn{iOv}, sub_fn{jOv}, mat2str(sv));
                    else
                        fprintf('  config.%s.%s = <%s>\n', fn{iOv}, sub_fn{jOv}, class(sv));
                    end
                end
            end
        else
            config.(fn{iOv}) = val;
            if ischar(val) || isstring(val)
                fprintf('  config.%s = %s\n', fn{iOv}, char(val));
            elseif isnumeric(val) && isscalar(val)
                fprintf('  config.%s = %g\n', fn{iOv}, val);
            else
                fprintf('  config.%s = <%s>\n', fn{iOv}, class(val));
            end
        end
    end
end

%% ========================================================================
%  DISCOVER FILES AND GROUP BY SERIES
%  ========================================================================
addpath(genpath(config.vadaScriptsPath));

pipelineDir = fileparts(mfilename('fullpath'));
if isempty(pipelineDir), pipelineDir = pwd; end
addpath(pipelineDir);
if exist('process_single_block', 'file') ~= 2
    error(['Cannot find process_single_block.m. Make sure all LAT_ULM .m files ' ...
           'are in the same folder and cd there or addpath it.']);
end
if ~exist(config.outputFolder, 'dir'), mkdir(config.outputFolder); end

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
        
        if isempty(vadaFiles)
            error('No %s files found in %s', config.modeName, config.dataFolder);
        end
        
        % Extract base filenames and read metadata
        nFound = numel(vadaFiles);
        fileList    = cell(nFound, 1);
        fileFolders = cell(nFound, 1);
        studyNames  = cell(nFound, 1);
        seriesNames = cell(nFound, 1);
        fileSizesGB = zeros(nFound, 1);
        
        fprintf('Reading metadata from %d file(s)...\n', nFound);
        for i = 1:nFound
            fullName = vadaFiles(i).name;
            baseName = fullName(1:end-numel(config.modeName));
            fileList{i}    = baseName;
            fileFolders{i} = vadaFiles(i).folder;
            fileSizesGB(i) = vadaFiles(i).bytes / 1e9;
            
            % Read study/series from XML
            xmlPath = fullfile(fileFolders{i}, [baseName config.modeName '.xml']);
            params = read_vada_xml_params(xmlPath);
            studyNames{i}  = get_param(params, 'Study-Name', 'UnknownStudy');
            seriesNames{i} = get_param(params, 'Series-Name', 'UnknownSeries');
        end
        
        % Deduplicate (same base name from different extensions)
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

%% ========================================================================
%  SERIES SUMMARY AND SELECTION (folder mode only)
%  ========================================================================
if strcmpi(config.mode, 'folder')
    % Build series groups: "StudyName | SeriesName"
    seriesKeys = strcat(studyNames, ' | ', seriesNames);
    uniqueSeries = unique(seriesKeys, 'stable');
    nSeries = numel(uniqueSeries);
    
    % Estimate processing time per block (~minutes based on typical throughput)
    % Rough estimate: ~2-5 min per block on A6000 depending on grid size
    estMinPerBlock = 3;  % Adjust after first run
    
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
    
    % Interactive selection
    fprintf('Options:\n');
    fprintf('  all     - Process all series\n');
    fprintf('  1,3,5   - Process specific series by number\n');
    fprintf('  1-3     - Process a range of series\n');
    fprintf('  q       - Quit\n\n');
    
    userInput = input('Select series to process: ', 's');
    userInput = strtrim(userInput);
    
    if strcmpi(userInput, 'q')
        fprintf('Cancelled.\n'); return;
    elseif strcmpi(userInput, 'all')
        selectedSeries = 1:nSeries;
    else
        % Parse comma-separated and range notation
        selectedSeries = [];
        parts = strsplit(userInput, ',');
        for p = 1:numel(parts)
            tok = strtrim(parts{p});
            rangeParts = strsplit(tok, '-');
            if numel(rangeParts) == 2
                r1 = str2double(rangeParts{1});
                r2 = str2double(rangeParts{2});
                if ~isnan(r1) && ~isnan(r2)
                    selectedSeries = [selectedSeries, r1:r2]; %#ok
                end
            else
                val = str2double(tok);
                if ~isnan(val)
                    selectedSeries = [selectedSeries, val]; %#ok
                end
            end
        end
        selectedSeries = unique(selectedSeries);
        selectedSeries = selectedSeries(selectedSeries >= 1 & selectedSeries <= nSeries);
    end
    
    if isempty(selectedSeries)
        fprintf('No valid series selected. Exiting.\n'); return;
    end
    
    % Filter to selected series only
    selectedMask = false(nFound, 1);
    for s = selectedSeries
        selectedMask = selectedMask | strcmp(seriesKeys, uniqueSeries{s});
    end
    
    fileList    = fileList(selectedMask);
    fileFolders = fileFolders(selectedMask);
    studyNames  = studyNames(selectedMask);
    seriesNames = seriesNames(selectedMask);
    seriesKeys  = seriesKeys(selectedMask);
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
                if ~isnan(r1) && ~isnan(r2), selectedBlocks = [selectedBlocks, r1:r2]; end %#ok
            else
                val = str2double(tok);
                if ~isnan(val), selectedBlocks = [selectedBlocks, val]; end %#ok
            end
        end
        selectedBlocks = unique(selectedBlocks);
        selectedBlocks = selectedBlocks(selectedBlocks >= 1 & selectedBlocks <= numel(fileList));
        
        if isempty(selectedBlocks)
            fprintf('No valid blocks. Exiting.\n'); return;
        end
        
        fileList    = fileList(selectedBlocks);
        fileFolders = fileFolders(selectedBlocks);
        studyNames  = studyNames(selectedBlocks);
        seriesNames = seriesNames(selectedBlocks);
        seriesKeys  = seriesKeys(selectedBlocks);
        fprintf('\nFiltered to %d block(s).\n', numel(fileList));
    end
end

numBlocks = numel(fileList);
fprintf('\nWill process %d block(s).\n', numBlocks);

%% ========================================================================
%  GPU SETUP
%  ========================================================================
if config.useGPU
    if gpuDeviceCount == 0
        warning('No GPU found, falling back to CPU.'); config.useGPU = false;
    else
        g = gpuDevice; reset(g);
        fprintf('\nGPU: %s (%.1f GB VRAM)\n', g.Name, g.TotalMemory/1e9);
    end
end

fprintf('\n=== LAT-ULM Pipeline (GPU=%s, mode=%s, %d blocks) ===\n', ...
    string(config.useGPU), config.mode, numBlocks);
pipelineTimer = tic;

%% ========================================================================
%  STEP 0: LOAD METADATA FROM FIRST BLOCK
%  ========================================================================
fprintf('\n[Step 0] Loading metadata from first block...\n');

% --- Auto-detect event structure if eventsPerFrame not specified ---
if isempty(config.eventsPerFrame)
    numProbe = 30;  % Generous sample to detect compound frame pattern
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
    % Unique signature: angle * 10 + polarity
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

    % Determine if PI is present: check if each unique angle appears more than once
    firstFrameAngles = probeAngles(1:eventsPerFrame);
    firstFramePolar  = probePolar(1:eventsPerFrame);
    uAngles = unique(firstFrameAngles, 'stable');
    hasPI = true;
    for a = 1:numel(uAngles)
        nOccur = sum(firstFrameAngles == uAngles(a));
        if nOccur < 2
            hasPI = false;
            break;
        end
    end
    config.hasPI = hasPI;
    config.numAngles = numel(uAngles);

    fprintf('  Auto-detected: %d events/frame, %d angles, PI=%s\n', ...
        eventsPerFrame, config.numAngles, string(hasPI));

    % Keep only the first eventsPerFrame events for metadata extraction
    VadaTest = VadaProbe(1:eventsPerFrame);
    clear VadaProbe;
else
    [VadaTest, Param, TxrParam, Config] = VsiVadaDataRead(...
        fileFolders{1}, fileList{1}, 1:config.eventsPerFrame, config.modeName);

    % Still need to detect PI and numAngles from the loaded events
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

% Pitch: VsiVadaDataRead documents ArrayPitch as micrometers, but for some
% probes (L38xp) VsiParseXml returns the value in mm (0.3, not 300).
% Detect unit by magnitude: real probe pitch is 0.05-1.0 mm (50-1000 um).
rawPitch = TxrParam.ArrayPitch;
fs_MHz = Param.SampleFreq;
depthOffset_mm = Param.DepthOffset;

if rawPitch == 0 || isnan(rawPitch)
    pitch_mm = config.pitchOverride;
    fprintf('  Pitch: %.4f mm (override, ArrayPitch=0)\n', pitch_mm);
elseif rawPitch < 10  % Already in mm
    pitch_mm = rawPitch;
    fprintf('  Pitch: %.4f mm (from metadata, interpreted as mm)\n', pitch_mm);
else  % In micrometers
    pitch_mm = rawPitch / 1000;
    fprintf('  Pitch: %.4f mm (converted from %.0f um)\n', pitch_mm, rawPitch);
end

% Speed of sound: XML "Vada-Mode/Speed-Of-Sound-Media" reflects the F2
% system setting (often 1480 m/s), not necessarily the actual medium.
% For agarose phantoms, true SoS ~ 1540 m/s. Override when needed.
c_xml = Param.SoSMedia;
if ~isempty(config.sosOverride)
    c = config.sosOverride;
    fprintf('  SoS: %.0f m/s (override; XML=%0.f m/s)\n', c, c_xml);
else
    c = c_xml;
    if c == 0, c = 1540; fprintf('  SoS: 1540 m/s (fallback, XML=0)\n');
    else, fprintf('  SoS: %.0f m/s (from XML)\n', c); end
end

fprintf('  %s | %d elem | pitch=%.4f mm (aperture=%.1f mm) | Fs=%.1f MHz\n', ...
    TxrParam.Name, TxrParam.ArrayNumElements, pitch_mm, ...
    (TxrParam.ArrayNumElements-1)*pitch_mm, fs_MHz);

% Read TX voltage from XML for blanking margin auto-scaling
if isempty(config.blankVoltage)
    xmlPath = fullfile(fileFolders{1}, [fileList{1} config.modeName '.xml']);
    xmlParams = read_vada_xml_params(xmlPath);
    voltHi = get_param(xmlParams, 'Vada-Mode/Voltage-Rail-High', '');
    voltLo = get_param(xmlParams, 'Vada-Mode/Voltage-Rail-Low', '');
    if ~isempty(voltHi)
        txVoltage = str2double(voltHi);
        fprintf('  TX voltage: %s-%s%% (from XML, using high=%s for blanking)\n', voltLo, voltHi, voltHi);
    else
        txVoltage = [];
        fprintf('  TX voltage: not found in XML (blanking uses base margin only)\n');
    end
else
    txVoltage = config.blankVoltage;
    fprintf('  TX voltage: %.0f%% (manual override)\n', txVoltage);
end
config.txVoltage = txVoltage;  % Store for downstream use

% Detect event structure
angleOrder = zeros(config.eventsPerFrame, 1);
polarityOrder = zeros(config.eventsPerFrame, 1);
for ev = 1:config.eventsPerFrame
    if isfield(VadaTest(ev).TxDelay, 'angle'), angleOrder(ev) = VadaTest(ev).TxDelay.angle; end
    if isfield(VadaTest(ev).Waveform, 'Channel') && isfield(VadaTest(ev).Waveform.Channel(1), 'invert')
        polarityOrder(ev) = VadaTest(ev).Waveform.Channel(1).invert;
    end
end
for ev = 1:config.eventsPerFrame
    fprintf('    Evt %d: angle=%+5.1f, inv=%d, Rx=[%d..%d](%d)\n', ev, ...
        angleOrder(ev), polarityOrder(ev), ...
        min(VadaTest(ev).Elements), max(VadaTest(ev).Elements), numel(VadaTest(ev).Elements));
end

uniqueAngles = unique(angleOrder, 'stable');
anglePairs = struct('angle',{},'posIdx',{},'negIdx',{},'rxElements',{});
for a = 1:numel(uniqueAngles)
    anglePairs(a).angle = uniqueAngles(a);
    idxs = find(angleOrder == uniqueAngles(a));
    if numel(idxs) >= 2
        % PI pair: two events per angle (positive + negative polarity)
        if polarityOrder(idxs(1)) == 0
            anglePairs(a).posIdx = idxs(1); anglePairs(a).negIdx = idxs(2);
        else
            anglePairs(a).posIdx = idxs(2); anglePairs(a).negIdx = idxs(1);
        end
    else
        % Single event per angle (no PI)
        anglePairs(a).posIdx = idxs(1);
        anglePairs(a).negIdx = [];
    end
    anglePairs(a).rxElements = VadaTest(anglePairs(a).posIdx).Elements;
end
config.numAngles = numel(uniqueAngles);

% Read TX frequency before clearing VadaTest
txFreq_MHz = 6;  % Default
if isfield(VadaTest(1).Waveform, 'Channel') && ~isempty(VadaTest(1).Waveform.Channel)
    if isfield(VadaTest(1).Waveform.Channel(1), 'frequency')
        txFreq_MHz = VadaTest(1).Waveform.Channel(1).frequency;
    end
end
lambda_mm = c * 1e-3 / txFreq_MHz;
nRxElements = numel(anglePairs(1).rxElements);
nSamplesPerEvent = size(VadaTest(1).Data, 1);

clear VadaTest;

%% ========================================================================
%  STEP 1: PRECOMPUTE DELAY TABLES (GPU)
%  ========================================================================
fprintf('\n[Step 1] Precomputing delay tables...\n');
elemPos_mm = ((1:TxrParam.ArrayNumElements) - (TxrParam.ArrayNumElements+1)/2) * pitch_mm;

% Auto-scale parameters from probe metadata if not manually set
rxElemPos = elemPos_mm(anglePairs(1).rxElements);
rxSpan = max(rxElemPos) - min(rxElemPos);
maxDataDepth_mm = depthOffset_mm + (nSamplesPerEvent / (fs_MHz * 2)) * (c * 1e-3);

fprintf('  TX freq: %.1f MHz, lambda: %.3f mm\n', txFreq_MHz, lambda_mm);
fprintf('  RX aperture: %.1f mm, depth range: %.1f-%.1f mm\n', rxSpan, depthOffset_mm, maxDataDepth_mm);

if isempty(config.bf.dx),     config.bf.dx = round(lambda_mm / 5, 4);          end
if isempty(config.bf.dz),     config.bf.dz = round(lambda_mm / 10, 4);         end
if isempty(config.bf.xRange)
    config.bf.xRange = [min(rxElemPos) - rxSpan*0.2, max(rxElemPos) + rxSpan*0.2];
end
if isempty(config.bf.zRange)
    config.bf.zRange = [max(depthOffset_mm, 0.5), maxDataDepth_mm];
end
if isempty(config.det.minSep_mm), config.det.minSep_mm = round(lambda_mm, 3); end
if isempty(config.sr.pixelSize_um), config.sr.pixelSize_um = max(1, round(lambda_mm / 5 * 1000)); end

fprintf('  Grid: dx=%.4f mm, dz=%.4f mm\n', config.bf.dx, config.bf.dz);
fprintf('  X range: [%.1f, %.1f] mm, Z range: [%.1f, %.1f] mm\n', ...
    config.bf.xRange(1), config.bf.xRange(2), config.bf.zRange(1), config.bf.zRange(2));
fprintf('  Detection minSep: %.3f mm, SR pixel: %d um\n', config.det.minSep_mm, config.sr.pixelSize_um);

xGrid = config.bf.xRange(1):config.bf.dx:config.bf.xRange(2);
zGrid = config.bf.zRange(1):config.bf.dz:config.bf.zRange(2);
nX = numel(xGrid); nZ = numel(zGrid);
fprintf('  Grid: %d x %d = %d pixels\n', nX, nZ, nX*nZ);

if config.zeroOnly
    % Find the 0-degree angle pair
    zeroAngleIdx = find([anglePairs.angle] == 0);
    if isempty(zeroAngleIdx)
        [~, zeroAngleIdx] = min(abs([anglePairs.angle]));
        fprintf('  WARNING: No exact 0deg angle. Using closest: %+.1f deg\n', anglePairs(zeroAngleIdx).angle);
    end
    config.zeroAngleIdx = zeroAngleIdx;
    config.zeroEventIdx = anglePairs(zeroAngleIdx).posIdx;  % Positive polarity event within each group
    
    % Only need one delay table
    rxPos_mm = elemPos_mm(anglePairs(zeroAngleIdx).rxElements);
    delayTables = {beamform_planewave_gpu([], rxPos_mm, anglePairs(zeroAngleIdx).angle, ...
        xGrid, zGrid, fs_MHz, c, depthOffset_mm, [])};
    fprintf('  Mode: ZERO-ONLY (0deg, no PI, complex IQ)\n');
    fprintf('  Using event %d from each group of %d\n', config.zeroEventIdx, config.eventsPerFrame);
    fprintf('  Effective frame rate: ~%.0f Hz\n', 1000 / (config.eventsPerFrame * 0.150));
else
    delayTables = cell(config.numAngles, 1);
    for a = 1:config.numAngles
        rxPos_mm = elemPos_mm(anglePairs(a).rxElements);
        delayTables{a} = beamform_planewave_gpu([], rxPos_mm, anglePairs(a).angle, ...
            xGrid, zGrid, fs_MHz, c, depthOffset_mm, []);
        fprintf('  Angle %+5.1f: ready\n', anglePairs(a).angle);
    end
    if config.hasPI
        fprintf('  Mode: ALL ANGLES with PI + compound\n');
    else
        fprintf('  Mode: ALL ANGLES, no PI (single event per angle)\n');
    end
    if config.blankSteering && config.hasPI
        pitch_for_blank = elemPos_mm(2) - elemPos_mm(1);
        if isempty(txVoltage)
            voltStr = 'unknown';
        else
            voltStr = sprintf('%.0f%%', txVoltage);
        end
        fprintf('  RF blanking enabled (base margin=%.1fx, voltage=%s):\n', ...
            config.blankMargin, voltStr);
        for a = 1:config.numAngles
            bi = compute_steering_blanking(anglePairs(a).angle, ...
                numel(anglePairs(a).rxElements), pitch_for_blank, c, fs_MHz, ...
                config.blankMargin, txVoltage);
            fprintf('    Angle %+5.1f: blank %d samples (%.1f us, margin=%.1fx, minDepth=%.2f mm)\n', ...
                bi.angle_deg, bi.nBlank, bi.delaySpread_us, bi.margin, bi.minDepth_mm);
        end
    end
end

%% ========================================================================
%  STEP 1b: COMPUTE BACKGROUND (if specified)
%  ========================================================================
bgMeanIQ = [];  % Empty = no background subtraction

if ~isempty(config.bgFile)
    fprintf('\n[Step 1b] Computing background from reference block...\n');
    bgFolder = config.bgFolder;
    if isempty(bgFolder), bgFolder = config.dataFolder; end

    fprintf('  Background file: %s\n', config.bgFile);

    try
        % Load background block metadata
        [VadaBg, ~, ~, BgConfig] = VsiVadaDataRead(bgFolder, config.bgFile, ...
            1:config.eventsPerFrame, config.modeName);
        numBgEvents = numel(BgConfig.PulseSequences(1).Events);
        numBgFrames = floor(numBgEvents / config.eventsPerFrame);
        clear VadaBg;

        maxBgFrames = min(numBgFrames, 500);
        numBgEventsLoad = maxBgFrames * config.eventsPerFrame;

        fprintf('  Loading %d events (%d frames)...\n', numBgEventsLoad, maxBgFrames);
        [VadaBg, ~, ~, ~] = VsiVadaDataRead(bgFolder, config.bgFile, ...
            1:numBgEventsLoad, config.modeName);

        % Beamform background into full IQ stack (needed for motion correction)
        fprintf('  Beamforming background...\n');
        bgStack = zeros(nZ, nX, maxBgFrames, 'single');

        wb = waitbar(0, 'Computing background...', 'Name', 'Background');
        wbClean = onCleanup(@() close_if_valid(wb));

        % Select beamformer for background computation
        if isfield(config, 'useCUDA') && config.useCUDA
            bg_bf_fn = @beamform_cuda;
        else
            bg_bf_fn = @beamform_planewave_gpu;
        end

        for iFrame = 1:maxBgFrames
            if mod(iFrame, 100) == 0
                waitbar(iFrame/maxBgFrames, wb, sprintf('BG frame %d/%d', iFrame, maxBgFrames));
            end
            baseEvt = (iFrame-1) * config.eventsPerFrame;

            if config.zeroOnly
                evIdx = baseEvt + config.zeroEventIdx;
                rfData = single(VadaBg(evIdx).Data);
                rxPos_mm = elemPos_mm(anglePairs(config.zeroAngleIdx).rxElements);
                bfImg = bg_bf_fn(rfData, rxPos_mm, ...
                    anglePairs(config.zeroAngleIdx).angle, ...
                    xGrid, zGrid, fs_MHz, c, depthOffset_mm, delayTables{1});
            else
                bfImg = complex(zeros(nZ, nX, 'single'));
                for a = 1:config.numAngles
                    if config.hasPI
                        rfData = single(VadaBg(baseEvt + anglePairs(a).posIdx).Data) + ...
                                 single(VadaBg(baseEvt + anglePairs(a).negIdx).Data);
                    else
                        rfData = single(VadaBg(baseEvt + anglePairs(a).posIdx).Data);
                    end
                    rxPos_mm = elemPos_mm(anglePairs(a).rxElements);
                    bfA = bg_bf_fn(rfData, rxPos_mm, anglePairs(a).angle, ...
                        xGrid, zGrid, fs_MHz, c, depthOffset_mm, delayTables{a});
                    bfImg = bfImg + single(bfA);
                end
            end
            bgStack(:,:,iFrame) = single(bfImg);
        end
        close_if_valid(wb);
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
        fprintf('  Background computed: mean of %d frames\n', maxBgFrames);
        fprintf('  Background signal level: %.1f (max envelope)\n', max(abs(bgMeanIQ(:))));

    catch ME
        fprintf('  WARNING: Background computation failed: %s\n', ME.message);
        fprintf('  Proceeding without background subtraction.\n');
        bgMeanIQ = [];
    end
else
    fprintf('\n[Step 1b] No background file specified (config.bgFile is empty)\n');
end

%% ========================================================================
%  STEP 2: BLOCK QUALITY ASSESSMENT (auto-discriminate bubble vs no-bubble)
%  ========================================================================
fprintf('\n[Step 2] Assessing block quality (screening for bubble presence)...\n');
fprintf('  Sampling %d frames per block for quick assessment...\n\n', 50);

blockQC_cells = cell(numBlocks, 1);
wb = waitbar(0, 'Assessing blocks...', 'Name', 'Block Quality');
wbCleanup0 = onCleanup(@() close_if_valid(wb));

for i = 1:numBlocks
    waitbar(i/numBlocks, wb, sprintf('Assessing block %d/%d: %s', i, numBlocks, fileList{i}));
    blockQC_cells{i} = assess_block_quality(fileFolders{i}, fileList{i}, config.modeName, ...
        config, anglePairs, elemPos_mm, delayTables, xGrid, zGrid, nX, nZ, ...
        fs_MHz, c, depthOffset_mm, bgMeanIQ);
end
close_if_valid(wb);
blockQC = [blockQC_cells{:}];
close_if_valid(wb);

% Print quality summary
fprintf('  =========================================================================\n');
fprintf('  BLOCK QUALITY ASSESSMENT\n');
fprintf('  =========================================================================\n');
fprintf('  %-4s %-35s %8s %8s %8s %8s %6s\n', ...
    '#', 'Filename', 'StdCntr', 'SNRdrop', 'Locs@5', 'Locs@8', 'Bubble');
fprintf('  %s\n', repmat('-', 1, 85));

for i = 1:numBlocks
    q = blockQC(i);
    bubbleStr = 'NO';
    if q.hasBubbles, bubbleStr = '>>> YES'; end
    fprintf('  %-4d %-35s %8.1f %7.1fdB %8.1f %8.2f %s\n', ...
        i, q.filename, q.stdContrast, q.snrDrop, q.locsMild, q.locsStrict, bubbleStr);
end

nBubbleBlocks = sum([blockQC.hasBubbles]);
fprintf('  =========================================================================\n');
fprintf('  Blocks with bubbles: %d / %d\n', nBubbleBlocks, numBlocks);
fprintf('  =========================================================================\n');
if nBubbleBlocks == 0
    fprintf('  WARNING: No blocks passed bubble detection.\n');
    fprintf('  Processing all blocks anyway (empty blocks contribute nothing).\n');
end
fprintf('\n');

%% ========================================================================
%  STEPS 3-6: PROCESS BLOCKS GROUPED BY SERIES
%  ========================================================================

% Determine series grouping for selected blocks
if strcmpi(config.mode, 'folder')
    processSeriesKeys = unique(seriesKeys, 'stable');
else
    processSeriesKeys = {'single_file'};
    seriesKeys = {'single_file'};
end
nProcessSeries = numel(processSeriesKeys);

fprintf('\n[Steps 3-6] Processing %d block(s) across %d series...\n', numBlocks, nProcessSeries);

% Grand accumulation across all series
allLocalizations   = [];
allTracks          = {};
allTimestamps      = [];
allBlockResults    = cell(numBlocks, 1);
globalFrameOffset  = 0;
totalFrames        = 0;
blockCounter       = 0;

% Per-series results
seriesResults = struct('key', {}, 'localizations', {}, 'tracks', {}, ...
    'timestamps', {}, 'numFrames', {}, 'numBlocks', {});

wb = waitbar(0, 'Processing...', 'Name', 'LAT-ULM Pipeline');
wbCleanup = onCleanup(@() close_if_valid(wb));

for iSeries = 1:nProcessSeries
    seriesKey = processSeriesKeys{iSeries};
    seriesMembers = find(strcmp(seriesKeys, seriesKey));
    nBlocksThisSeries = numel(seriesMembers);
    
    fprintf('\n  ====== Series %d/%d: %s (%d blocks) ======\n', ...
        iSeries, nProcessSeries, seriesKey, nBlocksThisSeries);
    
    seriesLocs = [];
    seriesTracks = {};
    seriesTimestamps = [];
    seriesFrameOffset = 0;
    
    for iBlockInSeries = 1:nBlocksThisSeries
        iBlock = seriesMembers(iBlockInSeries);
        blockCounter = blockCounter + 1;
        
        waitbar((blockCounter-1)/numBlocks, wb, sprintf(...
            'Series %d/%d | Block %d/%d: %s | %d locs total', ...
            iSeries, nProcessSeries, iBlockInSeries, nBlocksThisSeries, ...
            fileList{iBlock}, size(allLocalizations,1)));
        
        fprintf('\n    --- Block %d/%d (overall %d/%d): %s ---\n', ...
            iBlockInSeries, nBlocksThisSeries, blockCounter, numBlocks, fileList{iBlock});
        
        blockResult = process_single_block(fileFolders{iBlock}, fileList{iBlock}, config.modeName, ...
            config, anglePairs, elemPos_mm, delayTables, xGrid, zGrid, nX, nZ, fs_MHz, c, depthOffset_mm, bgMeanIQ);
        
        allBlockResults{blockCounter} = blockResult;
        
        % Offset frame indices (globally unique)
        if ~isempty(blockResult.localizations)
            blockResult.localizations(:,4) = blockResult.localizations(:,4) + globalFrameOffset;
            allLocalizations = [allLocalizations; blockResult.localizations]; %#ok
            seriesLocs = [seriesLocs; blockResult.localizations]; %#ok
        end
        
        for iT = 1:numel(blockResult.tracks)
            t = blockResult.tracks{iT};
            t(:,4) = t(:,4) + globalFrameOffset;
            allTracks{end+1} = t; %#ok
            seriesTracks{end+1} = t; %#ok
        end
        
        if ~isempty(blockResult.timestamps)
            allTimestamps = [allTimestamps; blockResult.timestamps]; %#ok
            seriesTimestamps = [seriesTimestamps; blockResult.timestamps]; %#ok
        end
        
        globalFrameOffset = globalFrameOffset + blockResult.numFrames;
        totalFrames = totalFrames + blockResult.numFrames;
        seriesFrameOffset = seriesFrameOffset + blockResult.numFrames;
        
        % Per-block diagnostics
        if config.perBlockDiagnostics
            estFrameRate = 1000 / (config.eventsPerFrame * 0.150);  % Approximate
            save_block_diagnostics(config.outputFolder, fileList{iBlock}, blockCounter, ...
                blockResult, config, xGrid, zGrid, fs_MHz, c, estFrameRate, 'latulm');
        end
    end
    
    % Store per-series results
    seriesResults(iSeries).key = seriesKey;
    seriesResults(iSeries).localizations = seriesLocs;
    seriesResults(iSeries).tracks = seriesTracks;
    seriesResults(iSeries).timestamps = seriesTimestamps;
    seriesResults(iSeries).numFrames = seriesFrameOffset;
    seriesResults(iSeries).numBlocks = nBlocksThisSeries;
    
    fprintf('\n    Series "%s" complete: %d locs, %d tracks, %d frames\n', ...
        seriesKey, size(seriesLocs,1), numel(seriesTracks), seriesFrameOffset);
end

waitbar(1, wb, sprintf('Done: %d locs, %d tracks', ...
    size(allLocalizations,1), numel(allTracks)));
pause(0.5);
close_if_valid(wb);

fprintf('\n=== Processing complete ===\n');
fprintf('  Series processed: %d\n', nProcessSeries);
fprintf('  Total blocks: %d\n', numBlocks);
fprintf('  Total frames: %d\n', totalFrames);
fprintf('  Total localizations: %d\n', size(allLocalizations, 1));
fprintf('  Total tracks: %d\n', numel(allTracks));
if ~isempty(allTracks)
    trackLengths = cellfun(@(t) size(t,1), allTracks);
    fprintf('  Track lengths: min=%d, med=%d, max=%d, mean=%.1f\n', ...
        min(trackLengths), round(median(trackLengths)), max(trackLengths), mean(trackLengths));
end

% Summary tables
fprintf('\n  Series Summary:\n');
fprintf('  %-4s  %-45s  %6s  %8s  %8s\n', '#', 'Series', 'Blocks', 'Locs', 'Tracks');
fprintf('  %s\n', repmat('-', 1, 75));
for s = 1:nProcessSeries
    fprintf('  %-4d  %-45s  %6d  %8d  %8d\n', s, seriesResults(s).key, ...
        seriesResults(s).numBlocks, size(seriesResults(s).localizations,1), ...
        numel(seriesResults(s).tracks));
end

fprintf('\n  Block Summary:\n');
fprintf('  %-4s %-40s %8s %8s %8s %8s\n', '#', 'Filename', 'Frames', 'Locs', 'Tracks', 'Time(s)');
fprintf('  %s\n', repmat('-', 1, 80));
for i = 1:blockCounter
    br = allBlockResults{i};
    fprintf('  %-4d %-40s %8d %8d %8d %8.1f\n', ...
        i, br.filename, br.numFrames, br.numLocs, br.numTracks, br.processingTime);
end

%% ========================================================================
%  STEP 7: SUPER-RESOLUTION RENDERING
%  ========================================================================
fprintf('\n[Step 7] Rendering super-resolution maps...\n');

srPx = config.sr.pixelSize_um / 1000;
srX = config.bf.xRange(1):srPx:config.bf.xRange(2);
srZ = config.bf.zRange(1):srPx:config.bf.zRange(2);
nSrX = numel(srX); nSrZ = numel(srZ);

densityMap = zeros(nSrZ, nSrX, 'single');
velMapX = zeros(nSrZ, nSrX, 'single');
velMapZ = zeros(nSrZ, nSrX, 'single');
velCount = zeros(nSrZ, nSrX, 'single');

if numel(allTimestamps) > 1
    frameRate = 1 / median(diff(allTimestamps) / 1000);
else
    frameRate = 1111;
end
fprintf('  Frame rate: %.1f Hz | SR pixel: %d um | Grid: %dx%d\n', ...
    frameRate, config.sr.pixelSize_um, nSrX, nSrZ);

nTracks = numel(allTracks);
if nTracks > 0
    wb2 = waitbar(0, 'Rendering...', 'Name', 'SR Rendering');
    wb2Cleanup = onCleanup(@() close_if_valid(wb2));
    
    for iT = 1:nTracks
        if mod(iT, 500) == 0 && ishandle(wb2)
            waitbar(iT/nTracks, wb2, sprintf('Tracks: %d/%d', iT, nTracks));
        end
        t = allTracks{iT};
        for p = 1:size(t, 1)
            [~,xi] = min(abs(srX - t(p,1)));
            [~,zi] = min(abs(srZ - t(p,2)));
            densityMap(zi,xi) = densityMap(zi,xi) + 1;
            if p > 1
                dt = (t(p,4) - t(p-1,4)) / frameRate;
                if dt > 0
                    velMapX(zi,xi) = velMapX(zi,xi) + (t(p,1)-t(p-1,1))/dt;
                    velMapZ(zi,xi) = velMapZ(zi,xi) + (t(p,2)-t(p-1,2))/dt;
                    velCount(zi,xi) = velCount(zi,xi) + 1;
                end
            end
        end
    end
    close_if_valid(wb2);
end

speedMap = zeros(size(densityMap), 'single');
m = velCount > 0;
speedMap(m) = sqrt((velMapX(m)./velCount(m)).^2 + (velMapZ(m)./velCount(m)).^2);

%% ========================================================================
%  STEP 8: FIGURES
%  ========================================================================
fprintf('\n[Step 8] Figures...\n');

figure('Name','Density','Position',[50 100 800 600]);
imagesc(srX, srZ, log10(densityMap+1)); axis image; colormap hot;
cb = colorbar; cb.Label.String = 'log_{10}(count+1)';
xlabel('Lateral [mm]'); ylabel('Axial [mm]');
title(sprintf('Density (%d blocks, %d tracks, %d locs, %d\\mum SR)', ...
    numBlocks, numel(allTracks), size(allLocalizations,1), config.sr.pixelSize_um));
saveas(gcf, fullfile(config.outputFolder, 'density_map.png'));
saveas(gcf, fullfile(config.outputFolder, 'density_map.fig'));

figure('Name','Velocity','Position',[50 100 800 600]);
imagesc(srX, srZ, speedMap); axis image; colormap jet;
cb = colorbar; cb.Label.String = 'Speed [mm/s]';
if any(speedMap(:)>0), caxis([0 prctile(speedMap(speedMap>0),95)]); end
xlabel('Lateral [mm]'); ylabel('Axial [mm]'); title('Velocity [mm/s]');
saveas(gcf, fullfile(config.outputFolder, 'velocity_map.png'));
saveas(gcf, fullfile(config.outputFolder, 'velocity_map.fig'));

figure('Name','Tracks','Position',[50 100 800 600]); hold on;
nPlot = min(numel(allTracks), 2000); 
if nPlot > 0
    cmap = jet(nPlot);
    idx = randperm(numel(allTracks), nPlot);
    for ii = 1:nPlot
        t = allTracks{idx(ii)};
        plot(t(:,1), t(:,2), '-', 'Color', [cmap(ii,:) 0.5], 'LineWidth', 0.5);
    end
end
set(gca,'YDir','reverse'); axis equal tight;
xlabel('Lateral [mm]'); ylabel('Axial [mm]');
title(sprintf('Tracks (%d/%d shown)', nPlot, numel(allTracks)));
saveas(gcf, fullfile(config.outputFolder, 'tracks.png'));

figure('Name','Stats','Position',[50 50 1400 400]);
subplot(1,4,1);
if ~isempty(allLocalizations), histogram(allLocalizations(:,3),50,'FaceColor',[0.2 0.4 0.8]); end
xlabel('Amplitude'); title('Detection Amplitudes');

subplot(1,4,2);
if ~isempty(allLocalizations)
    lpf = accumarray(allLocalizations(:,4), 1, [totalFrames 1]);
    plot(lpf, 'Color', [0.2 0.6 0.3]); xlabel('Frame'); title('Locs/Frame');
    fprintf('  Avg locs/frame: %.1f\n', mean(lpf));
end

subplot(1,4,3);
if ~isempty(allTracks)
    histogram(trackLengths, 50, 'FaceColor', [0.8 0.3 0.2]);
    xlabel('Frames'); title('Track Length');
end

subplot(1,4,4);
if ~isempty(allTracks)
    sp = [];
    for i = 1:numel(allTracks)
        t = allTracks{i};
        if size(t,1) > 1
            sp = [sp; sqrt(diff(t(:,1)).^2 + diff(t(:,2)).^2) * frameRate]; %#ok
        end
    end
    if ~isempty(sp)
        histogram(sp, 50, 'FaceColor', [0.6 0.2 0.6]);
        xlabel('mm/s'); title('Speed Distribution');
        fprintf('  Median speed: %.2f mm/s\n', median(sp));
    end
end
saveas(gcf, fullfile(config.outputFolder, 'statistics.png'));

%% ========================================================================
%  SAVE
%  ========================================================================
fprintf('\n[Saving]...\n');

R.config         = config;
R.fileList       = fileList;
R.fileFolders    = fileFolders;
R.blockQC        = blockQC;  % Quality assessment for all blocks
R.seriesResults  = seriesResults;
R.blockResults   = allBlockResults(1:blockCounter);
R.localizations  = allLocalizations;
R.tracks         = allTracks;
R.timestamps     = allTimestamps;
R.densityMap     = densityMap;
R.speedMap       = speedMap;
R.srX            = srX;
R.srZ            = srZ;
R.Param          = Param;
R.TxrParam       = TxrParam;
R.frameRate      = frameRate;
R.totalFrames    = totalFrames;
R.numBlocks      = numBlocks;
R.processingTime = toc(pipelineTimer);

save(fullfile(config.outputFolder, 'LAT_ULM_results.mat'), 'R', '-v7.3');

fprintf('\nSaved to %s\n', config.outputFolder);
fprintf('=== Done: %d series, %d blocks, %d frames, %d locs, %d tracks in %.1f min ===\n', ...
    nProcessSeries, numBlocks, totalFrames, size(allLocalizations,1), numel(allTracks), toc(pipelineTimer)/60);

if config.useGPU, reset(gpuDevice); end

%% ========================================================================
%  LOCAL HELPER FUNCTIONS
%  ========================================================================

function close_if_valid(h)
    if isvalid(h), close(h); end
end

function params = read_vada_xml_params(xmlPath)
% Read all <parameter name="..." value="..."/> from a .vada.xml file.
% Returns struct with field names derived from parameter names
% (hyphens/slashes replaced with underscores).
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
% Get parameter by original hyphenated name (auto-converts to struct field name).
fn = strrep(strrep(name, '-', '_'), '/', '_');
if isstruct(params) && isfield(params, fn)
    val = params.(fn);
else
    val = default;
end
end
