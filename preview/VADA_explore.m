%% VADA_EXPLORE.m
% Diagnostic script to inspect VADA data before running full LAT-ULM pipeline.
% GPU-accelerated beamforming for quick preview.
%
% SEQUENCE/TRANSDUCER AGNOSTIC: Auto-detects event structure (PI pairs,
% single events, any number of angles, any transducer).
%
% Run this FIRST to confirm:
%   1. Event structure (angles, polarities, apertures)
%   2. Pulse inversion quality (when PI is present)
%   3. Beamformed image preview
%   4. Noise floor and detection threshold estimates
%
% Author: Eli Wirth-Apley / Sun Lab, Northeastern University, February 2026

clearvars; close all; clc;

%% Configuration
dataFolder      = 'C:\path\to\VADA_data';
baseFilename    = '';  % CHANGE (no extension)
modeName        = '.vada';
vadaScriptsPath = 'C:\path\to\VevoF2\RFExampleScriptsF2\Matlab';
useGPU          = true;
sosOverride     = 1540;    % [m/s] Set to [] to use XML value, 1540 for agarose

% Angle mode (matches pipeline config.zeroOnly):
% true  = recommended: 0deg only, complex IQ, SVD-based clutter filter
% false = all angles with compound (PI if present)
zeroOnly        = true;

% Background subtraction: specify a water/PBS-only VADA block
% Leave empty to skip background subtraction
bgFile          = '';  % e.g., 'water_baseline_block' (no extension)
bgFolder        = '';  % Leave '' to use dataFolder

% Grid and detection parameters
% Set to [] for auto-scaling based on probe frequency, or override manually.
gridDx          = [];   % Lateral pixel [mm]. Auto: ~lambda/5
gridDz          = [];   % Axial pixel [mm].  Auto: ~lambda/10
gridXRange      = [];   % [min max] lateral [mm]. Auto from RX aperture
gridZRange      = [];   % [min max] axial [mm].   Auto from depth range
detMinSep       = [];   % Detection min separation [mm]. Auto: ~lambda
detROISize      = [];   % Detection ROI [px]. Auto: 7

addpath(genpath(vadaScriptsPath));

%% GPU check
if useGPU
    if gpuDeviceCount == 0
        fprintf('No GPU detected, using CPU.\n'); useGPU = false;
    else
        g = gpuDevice; reset(g);
        fprintf('GPU: %s (%.1f GB VRAM)\n', g.Name, g.TotalMemory/1e9);
    end
end

%% Step 1: Auto-detect event structure and load data
fprintf('\n[1/7] Auto-detecting event structure...\n');

% Load a generous sample to detect the compound frame structure
numProbe = 30;
[VadaProbe, Param, TxrParam, Config] = VsiVadaDataRead(...
    dataFolder, baseFilename, 1:numProbe, modeName);

% Extract angle and polarity from each event
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

% Create unique signature: angle * 10 + polarity
probeSig = probeAngles * 10 + probePolar;

% Find the period: when the full signature pattern first repeats
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

% Determine angles and whether PI is present
firstFrameAngles = probeAngles(1:eventsPerFrame);
firstFramePolar  = probePolar(1:eventsPerFrame);
uniqueAngles     = unique(firstFrameAngles, 'stable');
numAngles        = numel(uniqueAngles);

% Check if each angle appears more than once (PI pairs)
hasPI = true;
for a = 1:numAngles
    if sum(firstFrameAngles == uniqueAngles(a)) < 2
        hasPI = false;
        break;
    end
end

if hasPI
    numPairs = numAngles;  % Each angle has a PI pair
else
    numPairs = 0;  % No PI pairs
end

fprintf('  Detected: %d events/frame, %d unique angles, PI=%s\n', ...
    eventsPerFrame, numAngles, string(hasPI));

% Build angle info struct (works for PI and non-PI)
angleInfo = struct('angle',{},'posIdx',{},'negIdx',{},'rxElements',{});
for a = 1:numAngles
    angleInfo(a).angle = uniqueAngles(a);
    idxs = find(firstFrameAngles == uniqueAngles(a));
    if hasPI && numel(idxs) >= 2
        if firstFramePolar(idxs(1)) == 0
            angleInfo(a).posIdx = idxs(1);
            angleInfo(a).negIdx = idxs(2);
        else
            angleInfo(a).posIdx = idxs(2);
            angleInfo(a).negIdx = idxs(1);
        end
    else
        angleInfo(a).posIdx = idxs(1);
        angleInfo(a).negIdx = [];
    end
end

numLoadFrames = 3;
numLoad = numLoadFrames * eventsPerFrame;
fprintf('  Loading %d compound frames (%d events)...\n', numLoadFrames, numLoad);

clear VadaProbe;
[VadaMode, Param, TxrParam, Config] = VsiVadaDataRead(...
    dataFolder, baseFilename, 1:numLoad, modeName);

% Set rxElements from data
for a = 1:numAngles
    angleInfo(a).rxElements = VadaMode(angleInfo(a).posIdx).Elements;
end

%% Study/Series info from XML
fprintf('\n--- STUDY / SERIES ---\n');
xmlPath = fullfile(dataFolder, [baseFilename modeName '.xml']);
txVoltage = [];
if exist(xmlPath, 'file')
    xmlText = fileread(xmlPath);
    tokens = regexp(xmlText, '<parameter\s+name="([^"]+)"\s+value="([^"]*)"', 'tokens');
    xmlParams = struct();
    for ti = 1:numel(tokens)
        fn = strrep(strrep(tokens{ti}{1}, '-', '_'), '/', '_');
        xmlParams.(fn) = tokens{ti}{2};
    end
    if isfield(xmlParams, 'Study_Name'),  fprintf('  Study:  %s\n', xmlParams.Study_Name);  end
    if isfield(xmlParams, 'Series_Name'), fprintf('  Series: %s\n', xmlParams.Series_Name); end
    if isfield(xmlParams, 'Vada_Mode_User_Pulse_Sequence_Name')
        fprintf('  Sequence: %s\n', xmlParams.Vada_Mode_User_Pulse_Sequence_Name);
    end
    if isfield(xmlParams, 'Vada_Mode_Speed_Of_Sound_Media')
        fprintf('  SoS (XML): %s m/s\n', xmlParams.Vada_Mode_Speed_Of_Sound_Media);
    end
    if isfield(xmlParams, 'Vada_Mode_Voltage_Rail_Low') && isfield(xmlParams, 'Vada_Mode_Voltage_Rail_High')
        fprintf('  Voltage: %s-%s%%\n', xmlParams.Vada_Mode_Voltage_Rail_Low, xmlParams.Vada_Mode_Voltage_Rail_High);
        txVoltage = str2double(xmlParams.Vada_Mode_Voltage_Rail_High);
    end
    if isfield(xmlParams, 'Element_Pitch')
        fprintf('  Pitch (XML): %s mm\n', xmlParams.Element_Pitch);
    end
else
    fprintf('  WARNING: No .vada.xml found at %s\n', xmlPath);
end

%% System parameters
fprintf('\n--- SYSTEM PARAMETERS ---\n');
fprintf('  Transducer: %s\n', TxrParam.Name);
fprintf('  Total Elements: %d | Pitch (raw): %.4f | MaxTX: %d | MaxRX: %d\n', ...
    TxrParam.ArrayNumElements, TxrParam.ArrayPitch, TxrParam.MaxTxElements, TxrParam.MaxRxElements);
fprintf('  Fs: %.1f MHz | SoS: %.0f m/s | DepthOffset: %.2f mm\n', ...
    Param.SampleFreq, Param.SoSMedia, Param.DepthOffset);
fprintf('  RF per event: [%d samples x %d channels]\n', ...
    size(VadaMode(1).Data,1), size(VadaMode(1).Data,2));
rxElem = VadaMode(1).Elements;
nRx = numel(rxElem);
fprintf('  Active RX elements: [%d..%d] (%d elements)\n', ...
    min(rxElem), max(rxElem), nRx);

%% Step 2: Event structure
fprintf('\n[2/7] Event structure analysis...\n');
fprintf('  %-6s %-10s %-12s %-8s %-8s %-20s\n', ...
    'Event', 'Time(ms)', 'Angle', 'Invert', 'Freq', 'RxElements');

for ev = 1:min(2*eventsPerFrame, numel(VadaMode))
    angle_str = '?'; inv_str = '?'; freq_str = '?';

    if isfield(VadaMode(ev).TxDelay, 'angle')
        angle_str = sprintf('%+.1f deg', VadaMode(ev).TxDelay.angle);
    end
    if isfield(VadaMode(ev).Waveform, 'Channel')
        ch = VadaMode(ev).Waveform.Channel(1);
        if isfield(ch, 'invert'), inv_str = sprintf('%d', ch.invert); end
        if isfield(ch, 'frequency'), freq_str = sprintf('%.1fMHz', ch.frequency); end
    end

    fprintf('  %-6d %-10.3f %-12s %-8s %-8s [%d..%d](%d)\n', ...
        ev, VadaMode(ev).Timestamp, angle_str, inv_str, freq_str, ...
        min(VadaMode(ev).Elements), max(VadaMode(ev).Elements), numel(VadaMode(ev).Elements));
end

% Timing analysis
dt_events = diff([VadaMode(1:eventsPerFrame).Timestamp]);
fprintf('\n  Inter-event timing (ms): ');
fprintf('%.3f ', dt_events(1:min(5,end)));
if numel(VadaMode) > eventsPerFrame
    fprintf('\n  Frame period: %.3f ms -> %.0f Hz compound rate\n', ...
        VadaMode(eventsPerFrame+1).Timestamp - VadaMode(1).Timestamp, ...
        1000 / (VadaMode(eventsPerFrame+1).Timestamp - VadaMode(1).Timestamp));
end

%% Step 3: Pulse inversion quality (ONLY if PI is present)
if hasPI
    fprintf('\n[3/7] Pulse inversion quality check...\n');
    nPIRows = numPairs + 1;
    figH = min(300*nPIRows + 200, 1200);
    figure('Name', 'PI Quality', 'Position', [100 100 1400 figH]);
    midCh = ceil(size(VadaMode(1).Data, 2) / 2);

    for pair = 1:numPairs
        posIdx = angleInfo(pair).posIdx;
        negIdx = angleInfo(pair).negIdx;

        rfPos = double(VadaMode(posIdx).Data(:, midCh));
        rfNeg = double(VadaMode(negIdx).Data(:, midCh));
        rfSum = rfPos + rfNeg;

        corrVal = corrcoef(rfPos, rfNeg); corrVal = corrVal(1,2);
        cancDB = 10*log10(sum(rfSum.^2) / sum(rfPos.^2) + eps);

        % Time alignment check
        [xc, lags] = xcorr(rfPos, rfNeg, 10);
        [~, peakIdx] = max(abs(xc));
        lagAtPeak = lags(peakIdx);

        % Check if TX delay profiles match
        txDelayMatch = strcmp(VadaMode(posIdx).Event.txDelay, VadaMode(negIdx).Event.txDelay);

        % Per-element cancellation
        allCancDB = zeros(size(VadaMode(posIdx).Data, 2), 1);
        for ch = 1:size(VadaMode(posIdx).Data, 2)
            s = double(VadaMode(posIdx).Data(:,ch)) + double(VadaMode(negIdx).Data(:,ch));
            p = double(VadaMode(posIdx).Data(:,ch));
            allCancDB(ch) = 10*log10(sum(s.^2) / (sum(p.^2) + eps) + eps);
        end

        fprintf('  Pair %d (evt %d+%d, angle=%+.1f):\n', pair, posIdx, negIdx, ...
            angleInfo(pair).angle);
        fprintf('    Correlation: %.3f\n', corrVal);
        fprintf('    Cancellation: %.1f dB (mid ch), %.1f dB (mean all ch)\n', ...
            cancDB, mean(allCancDB));
        fprintf('    Cancellation range: [%.1f to %.1f] dB across elements\n', ...
            min(allCancDB), max(allCancDB));
        fprintf('    Time lag at peak xcorr: %d samples (%.1f ns)\n', ...
            lagAtPeak, lagAtPeak / (Param.SampleFreq) * 1000);
        fprintf('    TX delay profile match: %s\n', string(txDelayMatch));
        if cancDB < -20, fprintf('    Verdict: EXCELLENT\n');
        elseif cancDB < -10, fprintf('    Verdict: GOOD\n');
        else, fprintf('    Verdict: POOR - see diagnostics below\n'); end

        % Plots
        subplot(nPIRows, 3, (pair-1)*3 + 1);
        plot(rfPos,'b'); hold on; plot(rfNeg,'r','LineStyle','--');
        title(sprintf('Pair %d (%+.1f°): +/- RF', pair, angleInfo(pair).angle));
        legend('+','-'); xlabel('Sample'); ylabel('Amp');

        subplot(nPIRows, 3, (pair-1)*3 + 2);
        plot(rfSum,'k'); title(sprintf('PI Sum %.1fdB', cancDB));
        xlabel('Sample');

        subplot(nPIRows, 3, (pair-1)*3 + 3);
        plot(allCancDB, '.-'); hold on;
        yline(mean(allCancDB), 'r--', sprintf('mean=%.1f', mean(allCancDB)));
        xlabel('RX Channel'); ylabel('dB'); title('Per-element cancellation');
    end

    % Cross-correlation lag analysis
    if numPairs > 0
        subplot(nPIRows, 3, [numPairs*3+1, numPairs*3+2, numPairs*3+3]);
        colors = lines(numPairs);
        hold on;
        for pair = 1:numPairs
            posIdx = angleInfo(pair).posIdx;
            negIdx = angleInfo(pair).negIdx;
            nCh = size(VadaMode(posIdx).Data, 2);
            lagPerCh = zeros(nCh, 1);
            for ch = 1:nCh
                [xc, lags] = xcorr(double(VadaMode(posIdx).Data(:,ch)), ...
                                   double(VadaMode(negIdx).Data(:,ch)), 5);
                [~, pk] = max(abs(xc));
                lagPerCh(ch) = lags(pk);
            end
            plot(lagPerCh, '.-', 'Color', colors(pair,:), 'DisplayName', ...
                sprintf('%+.1f° (mean lag=%.2f)', angleInfo(pair).angle, mean(lagPerCh)));
        end
        legend('Location','best'); xlabel('RX Channel'); ylabel('Lag (samples)');
        title('PI Time Alignment: lag per channel per angle');
        yline(0, 'k--');
    end
    sgtitle('Pulse Inversion Diagnostics');
else
    fprintf('\n[3/7] No pulse inversion detected - skipping PI quality check.\n');
    fprintf('  Sequence has %d events/frame with %d unique angles (1 event per angle).\n', ...
        eventsPerFrame, numAngles);
end

%% Step 4: Beamforming diagnostics
fprintf('\n[4/7] Beamforming diagnostics (GPU=%s)...\n', string(useGPU));

% Pitch detection (probe-agnostic)
rawPitch = TxrParam.ArrayPitch;
if rawPitch == 0 || isnan(rawPitch)
    % Try XML first
    xmlPitchStr = '';
    if exist('xmlParams', 'var') && isfield(xmlParams, 'Element_Pitch')
        xmlPitchStr = xmlParams.Element_Pitch;
    end
    if ~isempty(xmlPitchStr)
        pitch_mm = str2double(xmlPitchStr);
        fprintf('  Pitch: %.4f mm (from XML Element-Pitch, ArrayPitch=0)\n', pitch_mm);
    else
        pitch_mm = 0.300;
        fprintf('  Pitch: %.3f mm (DEFAULT FALLBACK, ArrayPitch=0, no XML pitch)\n', pitch_mm);
        fprintf('  WARNING: Pitch was not found in metadata. Verify this value is correct.\n');
    end
elseif rawPitch < 10  % Already in mm
    pitch_mm = rawPitch;
    fprintf('  Pitch: %.4f mm (from metadata, interpreted as mm)\n', pitch_mm);
else  % In micrometers
    pitch_mm = rawPitch / 1000;
    fprintf('  Pitch: %.4f mm (converted from %.0f um)\n', pitch_mm, rawPitch);
end

elemPos_mm = ((1:TxrParam.ArrayNumElements) - (TxrParam.ArrayNumElements+1)/2) * pitch_mm;
rxPos_mm = elemPos_mm(rxElem);
fprintf('  Full array: %d elements, %.1f mm span\n', TxrParam.ArrayNumElements, ...
    (TxrParam.ArrayNumElements-1)*pitch_mm);
fprintf('  RX aperture: elements %d-%d (%d active), span=%.1f mm\n', ...
    min(rxElem), max(rxElem), nRx, max(rxPos_mm)-min(rxPos_mm));

c = Param.SoSMedia; if c==0, c=1540; end
if ~isempty(sosOverride)
    fprintf('  SoS: %.0f m/s (override; XML=%.0f m/s)\n', sosOverride, c);
    c = sosOverride;
end
fs_MHz = Param.SampleFreq;
depthOffset_mm = Param.DepthOffset;

% TX delay diagnostic
fprintf('\n  --- TX Delay Diagnostic ---\n');
for a = 1:numAngles
    posIdx = angleInfo(a).posIdx;
    txd = VadaMode(posIdx).TxDelay;
    angle_deg = txd.angle;

    fprintf('  Angle %+.1f: type=%s\n', angle_deg, txd.type);

    vadaDelays_ns = txd.delays;
    txElem = VadaMode(posIdx).Event.Tx.elements + VadaMode(posIdx).Event.Tx.offset;
    txElem = txElem(txElem > 0 & txElem <= TxrParam.ArrayNumElements);
    txPos_mm = elemPos_mm(txElem);
    idealDelays_ns = (txPos_mm * sind(angle_deg) / (c * 1e-6));
    idealDelays_ns = idealDelays_ns - min(idealDelays_ns);

    if numel(vadaDelays_ns) == numel(idealDelays_ns)
        vadaNorm = vadaDelays_ns - min(vadaDelays_ns);
        delayErr_ns = vadaNorm(:) - idealDelays_ns(:);
        fprintf('    Delay error: mean=%.1f ns, max=%.1f ns, std=%.1f ns\n', ...
            mean(delayErr_ns), max(abs(delayErr_ns)), std(delayErr_ns));
        fprintf('    (%.2f samples mean error at %.0f MHz)\n', ...
            mean(delayErr_ns) * fs_MHz / 1000, fs_MHz);
    else
        fprintf('    WARNING: VADA has %d delays, model expects %d\n', ...
            numel(vadaDelays_ns), numel(idealDelays_ns));
    end

    fprintf('    VADA delays (first 5): ');
    fprintf('%.1f ', vadaDelays_ns(1:min(5,end))); fprintf('ns\n');
end

% Auto-scale grid and detection parameters based on TX frequency
txFreq_MHz = 6;  % Default
if isfield(VadaMode(1).Waveform, 'Channel') && ~isempty(VadaMode(1).Waveform.Channel)
    if isfield(VadaMode(1).Waveform.Channel(1), 'frequency')
        txFreq_MHz = VadaMode(1).Waveform.Channel(1).frequency;
    end
end
lambda_mm = (c * 1e-3) / txFreq_MHz;
fprintf('\n  TX frequency: %.1f MHz, wavelength: %.3f mm\n', txFreq_MHz, lambda_mm);

if isempty(gridDx),     gridDx = round(lambda_mm / 5, 3);          end
if isempty(gridDz),     gridDz = round(lambda_mm / 10, 4);         end
if isempty(gridXRange)
    rxSpan = max(rxPos_mm) - min(rxPos_mm);
    gridXRange = [min(rxPos_mm) - rxSpan*0.2, max(rxPos_mm) + rxSpan*0.2];
    fprintf('  Lateral range from RX aperture: %.1f to %.1f mm\n', gridXRange(1), gridXRange(2));
end
if isempty(gridZRange)
    nSamples = size(VadaMode(1).Data, 1);
    maxDataDepth_mm = depthOffset_mm + (nSamples / (fs_MHz * 2)) * (c * 1e-3);
    gridZRange = [max(depthOffset_mm, 0.5), maxDataDepth_mm];
    fprintf('  Depth range from data: %.1f to %.1f mm (offset=%.1f, %d samples)\n', ...
        gridZRange(1), gridZRange(2), depthOffset_mm, nSamples);
end
if isempty(detMinSep),  detMinSep = round(lambda_mm, 3);            end
if isempty(detROISize), detROISize = 7;                              end

xGrid = gridXRange(1):gridDx:gridXRange(2);
zGrid = gridZRange(1):gridDz:gridZRange(2);

fprintf('  Grid: %d x %d pixels (dx=%.3f mm, dz=%.4f mm)\n', ...
    numel(xGrid), numel(zGrid), gridDx, gridDz);
fprintf('  Detection: minSep=%.3f mm, ROI=%d px\n', detMinSep, detROISize);

% Precompute delay tables for each angle
angles = [angleInfo.angle];
dtables = cell(numAngles, 1);
for a = 1:numAngles
    rxPosA = elemPos_mm(angleInfo(a).rxElements);
    dtables{a} = beamform_planewave_gpu([], rxPosA, angles(a), ...
        xGrid, zGrid, fs_MHz, c, depthOffset_mm, []);
end

% Beamform all events of first compound frame independently
bfEvents = zeros(numel(zGrid), numel(xGrid), eventsPerFrame);
wb = waitbar(0, 'Beamforming diagnostics...', 'Name', 'VADA Explore');
for ev = 1:eventsPerFrame
    waitbar(ev/eventsPerFrame, wb, sprintf('Beamforming event %d/%d...', ev, eventsPerFrame));
    rfData = single(VadaMode(ev).Data);
    % Find which angle this event belongs to
    evAngle = 0;
    if isfield(VadaMode(ev).TxDelay, 'angle')
        evAngle = VadaMode(ev).TxDelay.angle;
    end
    aIdx = find(angles == evAngle, 1);
    if isempty(aIdx), aIdx = 1; end
    rxPosEv = elemPos_mm(VadaMode(ev).Elements);
    bfEvents(:,:,ev) = beamform_planewave_gpu(rfData, rxPosEv, evAngle, ...
        xGrid, zGrid, fs_MHz, c, depthOffset_mm, dtables{aIdx});
end
close(wb);

%% === Compounding strategy comparison ===

% RF blanking (only meaningful with PI + steered angles)
if hasPI
    if isempty(txVoltage)
        fprintf('\n  RF blanking (base margin=1.5x, voltage=unknown):\n');
    else
        fprintf('\n  RF blanking (base margin=1.5x, voltage=%.0f%%):\n', txVoltage);
    end
    blankSamples = zeros(numAngles, 1);
    for a = 1:numAngles
        bi = compute_steering_blanking(angles(a), nRx, pitch_mm, c, fs_MHz, 1.5, txVoltage);
        blankSamples(a) = bi.nBlank;
        fprintf('    Angle %+5.1f: blank %d samples (%.1f us, eff margin=%.1fx)\n', ...
            angles(a), bi.nBlank, bi.delaySpread_us, bi.margin);
    end

    % Beamform blanked PI for each angle
    bfPI_blank = zeros(numel(zGrid), numel(xGrid), numAngles);
    bfBm_blank = zeros(numel(zGrid), numel(xGrid), numAngles);
    for a = 1:numAngles
        posIdx = angleInfo(a).posIdx;
        negIdx = angleInfo(a).negIdx;
        rfPos = single(VadaMode(posIdx).Data);
        rfNeg = single(VadaMode(negIdx).Data);
        nb = blankSamples(a);
        if nb > 0
            rfPos(1:nb,:) = 0;
            rfNeg(1:nb,:) = 0;
        end
        rxPosA = elemPos_mm(angleInfo(a).rxElements);
        bfPI_blank(:,:,a) = beamform_planewave_gpu(rfPos+rfNeg, rxPosA, angles(a), ...
            xGrid, zGrid, fs_MHz, c, depthOffset_mm, dtables{a});
        bfBm_blank(:,:,a) = beamform_planewave_gpu(rfPos-rfNeg, rxPosA, angles(a), ...
            xGrid, zGrid, fs_MHz, c, depthOffset_mm, dtables{a});
    end
end

% Find 0deg event index for display
zeroAngleIdx = find(angles == 0, 1);
if isempty(zeroAngleIdx), [~, zeroAngleIdx] = min(abs(angles)); end
zeroPosEvt = angleInfo(zeroAngleIdx).posIdx;

% Build comparison images
img_single = abs(bfEvents(:,:,zeroPosEvt));

% Incoherent compound: sum envelopes of all positive events
img_incoherent = zeros(numel(zGrid), numel(xGrid));
for a = 1:numAngles
    img_incoherent = img_incoherent + abs(bfEvents(:,:,angleInfo(a).posIdx));
end

% Coherent compound positive events only (no PI)
cohSum = complex(zeros(numel(zGrid), numel(xGrid)));
for a = 1:numAngles
    cohSum = cohSum + bfEvents(:,:,angleInfo(a).posIdx);
end
img_coh_pos = abs(cohSum);

% 0deg only (pipeline mode)
img_0deg_iq = abs(bfEvents(:,:,zeroPosEvt));

if hasPI
    % PI compound original
    piImgs = zeros(numel(zGrid), numel(xGrid), numAngles);
    for a = 1:numAngles
        piImgs(:,:,a) = bfEvents(:,:,angleInfo(a).posIdx) + bfEvents(:,:,angleInfo(a).negIdx);
    end
    img_piCompound = abs(sum(piImgs, 3));
    img_piCompound_blank = abs(sum(bfPI_blank, 3));

    % B-mode (PI diff) compound
    diffImgs = zeros(numel(zGrid), numel(xGrid), numAngles);
    for a = 1:numAngles
        diffImgs(:,:,a) = bfEvents(:,:,angleInfo(a).posIdx) - bfEvents(:,:,angleInfo(a).negIdx);
    end
    img_bmode = abs(sum(diffImgs, 3));
    img_bmode_blank = abs(sum(bfBm_blank, 3));

    % Plot 8-panel comparison (PI present)
    figure('Name', 'Compounding Comparison', 'Position', [50 50 1800 800]);
    titles = {'Single 0°', sprintf('Incoherent %d-angle', numAngles), ...
              sprintf('Coherent %d-pos (no PI)', numAngles), ...
              'PI compound (ORIGINAL)', 'PI compound (BLANKED)', ...
              'B-mode (ORIGINAL)', 'B-mode (BLANKED)', '0° IQ (PIPELINE)'};
    imgs = {img_single, img_incoherent, img_coh_pos, img_piCompound, ...
            img_piCompound_blank, img_bmode, img_bmode_blank, img_0deg_iq};

    for i = 1:8
        subplot(2, 4, i);
        img = imgs{i};
        imagesc(xGrid, zGrid, 20*log10(img / max(img(:)) + eps));
        axis image; colormap gray; colorbar; caxis([-60 0]);
        xlabel('Lat [mm]'); ylabel('Ax [mm]'); title(titles{i}, 'FontSize', 9);
    end
    sgtitle('Compounding: Original vs RF-Blanked PI');
else
    % Plot comparison without PI (fewer panels)
    figure('Name', 'Compounding Comparison', 'Position', [50 50 1400 500]);
    titles = {'Single 0°', sprintf('Incoherent %d-angle', numAngles), ...
              sprintf('Coherent %d-angle', numAngles), '0° IQ (PIPELINE)'};
    imgs = {img_single, img_incoherent, img_coh_pos, img_0deg_iq};

    for i = 1:numel(imgs)
        subplot(1, numel(imgs), i);
        img = imgs{i};
        imagesc(xGrid, zGrid, 20*log10(img / max(img(:)) + eps));
        axis image; colormap gray; colorbar; caxis([-60 0]);
        xlabel('Lat [mm]'); ylabel('Ax [mm]'); title(titles{i}, 'FontSize', 9);
    end
    sgtitle('Compounding Comparison (no PI)');
end

% SoS sensitivity (separate figure)
figure('Name', 'SoS Sensitivity', 'Position', [100 100 600 400]);
sosTest = [1480, 1510, 1540, 1570];
hold on;
midCol = round(numel(xGrid)/2);
testEvIdx = angleInfo(zeroAngleIdx).posIdx;
for si = 1:numel(sosTest)
    dt_test = beamform_planewave_gpu([], rxPos_mm, 0, ...
        xGrid, zGrid, fs_MHz, sosTest(si), depthOffset_mm, []);
    bf_test = beamform_planewave_gpu(single(VadaMode(testEvIdx).Data), rxPos_mm, 0, ...
        xGrid, zGrid, fs_MHz, sosTest(si), depthOffset_mm, dt_test);
    axLine = abs(bf_test(:, midCol));
    plot(zGrid, 20*log10(axLine / max(axLine) + eps), 'DisplayName', sprintf('%d m/s', sosTest(si)));
end
legend('Location', 'best'); xlabel('Depth [mm]'); ylabel('dB');
title('SoS sensitivity (0° axial)'); grid on;

%% Step 5: COMPLEX IQ SVD (0deg only - matches pipeline zeroOnly mode)
fprintf('\n[5/7] Complex IQ SVD (0deg, no PI)...\n');
fprintf('  This matches pipeline config.zeroOnly=true\n\n');

numCompoundTarget = 500;
numSvdEvents = min(numCompoundTarget * eventsPerFrame, numel(Config.PulseSequences(1).Events));
numCompound = floor(numSvdEvents / eventsPerFrame);
numSvdEvents = numCompound * eventsPerFrame;

% 0deg event position within each group
zeroEventInGroup = angleInfo(zeroAngleIdx).posIdx;
fprintf('  0deg event index within group: %d (of %d)\n', zeroEventInGroup, eventsPerFrame);

fprintf('  Loading %d events (%d compound frames)...\n', numSvdEvents, numCompound);

try
    [VadaSvd, ~, ~, ~] = VsiVadaDataRead(dataFolder, baseFilename, 1:numSvdEvents, modeName);
    fprintf('  Loaded. Beamforming 0deg events to complex IQ...\n');

    nZs = numel(zGrid); nXs = numel(xGrid);
    dt0 = dtables{zeroAngleIdx};

    iqStack = zeros(nZs, nXs, numCompound, 'single');

    wb = waitbar(0, 'Complex IQ beamforming (0deg)...', 'Name', 'SVD Explore');
    for iFrame = 1:numCompound
        if mod(iFrame, 50) == 0
            waitbar(iFrame/numCompound, wb, sprintf('Beamforming %d/%d...', iFrame, numCompound));
        end

        evIdx = (iFrame-1) * eventsPerFrame + zeroEventInGroup;
        rfData = single(VadaSvd(evIdx).Data);
        rxPosZ = elemPos_mm(angleInfo(zeroAngleIdx).rxElements);

        bfComplex = beamform_planewave_gpu(rfData, rxPosZ, angles(zeroAngleIdx), ...
            xGrid, zGrid, fs_MHz, c, depthOffset_mm, dt0);
        iqStack(:,:,iFrame) = single(bfComplex);
    end
    close(wb);
    clear VadaSvd;

    fprintf('  Complex IQ stack: [%d x %d x %d]\n', nZs, nXs, numCompound);

    % Background subtraction
    if ~isempty(bgFile)
        fprintf('  Computing background from: %s\n', bgFile);
        bgDir = bgFolder; if isempty(bgDir), bgDir = dataFolder; end

        try
            [VadaBg, ~, ~, BgCfg] = VsiVadaDataRead(bgDir, bgFile, ...
                1:min(500*eventsPerFrame, numel(Config.PulseSequences(1).Events)), modeName);
            nBgFrames = floor(numel(VadaBg) / eventsPerFrame);
            bgAccum = complex(zeros(nZs, nXs, 'single'));

            for bf = 1:nBgFrames
                bgEvIdx = (bf-1)*eventsPerFrame + zeroEventInGroup;
                bgRf = single(VadaBg(bgEvIdx).Data);
                rxPosBg = elemPos_mm(angleInfo(zeroAngleIdx).rxElements);
                bgBf = beamform_planewave_gpu(bgRf, rxPosBg, angles(zeroAngleIdx), ...
                    xGrid, zGrid, fs_MHz, c, depthOffset_mm, dt0);
                bgAccum = bgAccum + single(bgBf);
            end
            clear VadaBg;
            bgMeanIQ = bgAccum / nBgFrames;

            fprintf('  Background: mean of %d frames. Subtracting...\n', nBgFrames);
            for iFrame = 1:numCompound
                iqStack(:,:,iFrame) = iqStack(:,:,iFrame) - bgMeanIQ;
            end
            fprintf('  Background subtracted.\n');
        catch ME
            fprintf('  WARNING: Background failed: %s\n', ME.message);
        end
    end

    % SVD on COMPLEX Casorati matrix
    fprintf('  Computing SVD on complex data (%d frames)...\n', numCompound);
    Casorati = reshape(iqStack, nZs*nXs, numCompound);

    if useGPU
        [U, S, V] = svd(gpuArray(single(Casorati)), 'econ');
        U = gather(U); S = gather(S); V = gather(V);
    else
        [U, S, V] = svd(single(Casorati), 'econ');
    end
    singVals = diag(S);
    clear Casorati;

    fprintf('  SV ratio 1/2: %.1f, 1/5: %.1f, 1/10: %.1f, 1/20: %.1f\n', ...
        singVals(1)/singVals(2), singVals(1)/singVals(min(5,end)), ...
        singVals(1)/singVals(min(10,end)), singVals(1)/singVals(min(20,end)));

    % === Figure 1: SVD filtered mean images at different cutoffs ===
    svdCutoffs = [1, 2, 5, 10, 20, 50, 100];
    svdCutoffs = svdCutoffs(svdCutoffs < numCompound);

    figure('Name', 'Complex SVD Clutter Filter', 'Position', [50 50 1800 700]);
    nCols = ceil((numel(svdCutoffs)+2) / 2);

    % Mean envelope (tissue)
    subplot(2, nCols, 1);
    meanEnv = mean(abs(iqStack), 3);
    imagesc(xGrid, zGrid, 20*log10(meanEnv/max(meanEnv(:))+eps));
    axis image; colormap gray; colorbar; caxis([-60 0]);
    title('Mean envelope (tissue)'); xlabel('Lat [mm]'); ylabel('Ax [mm]');

    % Temporal std of COMPLEX data
    subplot(2, nCols, 2);
    complexStd = abs(std(iqStack, 0, 3));
    imagesc(xGrid, zGrid, 20*log10(complexStd/max(complexStd(:))+eps));
    axis image; colormap gray; colorbar; caxis([-60 0]);
    title('Complex temporal std'); xlabel('Lat [mm]'); ylabel('Ax [mm]');

    for ci = 1:numel(svdCutoffs)
        cutoff = svdCutoffs(ci);
        nKeep = size(U, 2);
        if cutoff < nKeep
            filtered = U(:, cutoff+1:nKeep) * S(cutoff+1:nKeep, cutoff+1:nKeep) * V(:, cutoff+1:nKeep)';
        else
            filtered = zeros(nZs*nXs, numCompound, 'single');
        end
        filtStack = reshape(filtered, nZs, nXs, numCompound);
        filtMeanEnv = mean(abs(filtStack), 3);

        subplot(2, nCols, ci+2);
        imagesc(xGrid, zGrid, 20*log10(filtMeanEnv/max(filtMeanEnv(:))+eps));
        axis image; colormap gray; colorbar; caxis([-60 0]);
        title(sprintf('SVD cut=%d', cutoff));
        xlabel('Lat [mm]'); ylabel('Ax [mm]');
    end
    sgtitle(sprintf('Complex IQ SVD (%d frames, 0deg, phase-preserved)', numCompound));

    % === Figure 2: Individual filtered frames ===
    figure('Name', 'Individual SVD-Filtered Frames', 'Position', [50 50 1600 800]);

    for testCut = [5, 20]
        if testCut >= size(U, 2), continue; end
        filtRecon = U(:, testCut+1:end) * S(testCut+1:end, testCut+1:end) * V(:, testCut+1:end)';
        filtFrames = reshape(filtRecon, nZs, nXs, numCompound);

        frameIdxs = round(linspace(1, numCompound, 4));
        rowOffset = (testCut == 20) * 4;
        for fi = 1:4
            subplot(2, 4, fi + rowOffset);
            frame = abs(filtFrames(:,:,frameIdxs(fi)));
            imagesc(xGrid, zGrid, frame);
            axis image; colormap hot; colorbar;
            title(sprintf('Fr %d (cut=%d)', frameIdxs(fi), testCut), 'FontSize', 9);
            xlabel('Lat [mm]'); ylabel('Ax [mm]');
        end
    end
    sgtitle('Individual SVD-filtered frames (top: cut=5, bottom: cut=20)');

    % === Figure 3: SV spectrum ===
    figure('Name', 'SVD Spectrum (Complex)', 'Position', [100 100 800 500]);

    subplot(1,2,1);
    semilogy(singVals / singVals(1), 'b.-', 'LineWidth', 1.5);
    hold on;
    for ci = 1:numel(svdCutoffs)
        xline(svdCutoffs(ci)+0.5, '--', sprintf('%d', svdCutoffs(ci)));
    end
    xlabel('SV Index'); ylabel('Normalized Magnitude');
    title(sprintf('Full spectrum (%d SVs)', numel(singVals)));
    grid on;

    subplot(1,2,2);
    nShow = min(50, numel(singVals));
    bar(singVals(1:nShow) / singVals(1), 'FaceColor', [0.3 0.5 0.8]);
    xlabel('SV Index'); ylabel('Normalized Magnitude');
    title('First 50 SVs (linear)');
    grid on;

    sgtitle('SVD Spectrum (Complex IQ, 0deg)');

    % SNR estimates
    fprintf('\n  SNR estimates by SVD cutoff:\n');
    for testCut = [2, 5, 10, 20]
        if testCut >= size(U, 2), continue; end
        filt = U(:, testCut+1:end) * S(testCut+1:end, testCut+1:end) * V(:, testCut+1:end)';
        fImg = reshape(filt, nZs, nXs, numCompound);
        fEnv = mean(abs(fImg), 3);
        nReg = fEnv(1:20, :);
        snr = 20*log10(max(fEnv(:)) / (std(nReg(:)) + eps));
        fprintf('    cut=%2d: SNR=%.1f dB, max=%.1f, noise_std=%.1f\n', ...
            testCut, snr, max(fEnv(:)), std(nReg(:)));
    end

    %% Step 5b: DETECTION TUNING
    fprintf('\n[5b/7] Detection parameter tuning...\n');

    tuneCutoff = 5;
    if tuneCutoff >= size(U, 2), tuneCutoff = min(2, size(U,2)-1); end

    filtRecon = U(:, tuneCutoff+1:end) * S(tuneCutoff+1:end, tuneCutoff+1:end) * V(:, tuneCutoff+1:end)';
    tuneStack = reshape(filtRecon, nZs, nXs, numCompound);

    dx = xGrid(2) - xGrid(1);
    dz = zGrid(2) - zGrid(1);

    testThresholds = [3, 5, 8, 10, 15, 20, 25, 30];
    locsPerFrame = zeros(numel(testThresholds), 1);

    detParams.method = 'threshold';
    detParams.minSep_mm = detMinSep;
    detParams.roiSize_px = detROISize;

    fprintf('  SVD cutoff=%d, minSep=%.1f mm\n', tuneCutoff, detParams.minSep_mm);
    fprintf('  %-12s %-15s %-15s\n', 'Threshold', 'Locs/frame', 'Total locs');
    fprintf('  %s\n', repmat('-', 1, 42));

    for ti = 1:numel(testThresholds)
        detParams.threshold = testThresholds(ti);
        totalLocs = 0;
        for f = 1:numCompound
            frame = abs(tuneStack(:,:,f));
            dets = detect_microbubbles(frame, detParams, dx, dz);
            totalLocs = totalLocs + size(dets, 1);
        end
        locsPerFrame(ti) = totalLocs / numCompound;
        fprintf('  %-12d %-15.1f %-15d\n', testThresholds(ti), locsPerFrame(ti), totalLocs);
    end

    % Plot threshold curve
    figure('Name', 'Detection Tuning', 'Position', [100 100 1400 500]);

    subplot(1,3,1);
    semilogy(testThresholds, locsPerFrame, 'b.-', 'LineWidth', 2, 'MarkerSize', 12);
    hold on;
    yline(3, 'g--', 'Target: 1-5/frame', 'LineWidth', 1.5);
    yline(1, 'g--');
    yline(5, 'g--');
    xlabel('Detection Threshold (noise std)');
    ylabel('Localizations per Frame');
    title('Threshold Tuning Curve');
    grid on;

    [~, bestIdx] = min(abs(locsPerFrame - 3));
    bestThresh = testThresholds(bestIdx);
    fprintf('\n  Recommended threshold for ~3 locs/frame: %d\n', bestThresh);

    detParams.threshold = bestThresh;
    sampleFrame = round(numCompound / 2);
    frame = abs(tuneStack(:,:,sampleFrame));
    dets = detect_microbubbles(frame, detParams, dx, dz);

    subplot(1,3,2);
    imagesc(xGrid, zGrid, frame);
    axis image; colormap hot; colorbar; hold on;
    if ~isempty(dets)
        detX = xGrid(1) + (dets(:,2)-1) * dx;
        detZ = zGrid(1) + (dets(:,1)-1) * dz;
        plot(detX, detZ, 'go', 'MarkerSize', 8, 'LineWidth', 2);
    end
    title(sprintf('Frame %d: %d detections (thresh=%d)', sampleFrame, size(dets,1), bestThresh));
    xlabel('Lat [mm]'); ylabel('Ax [mm]');

    quickDensity = zeros(nZs, nXs);
    for f = 1:numCompound
        frame = abs(tuneStack(:,:,f));
        dets = detect_microbubbles(frame, detParams, dx, dz);
        if ~isempty(dets)
            for d = 1:size(dets, 1)
                quickDensity(dets(d,1), dets(d,2)) = quickDensity(dets(d,1), dets(d,2)) + 1;
            end
        end
    end

    subplot(1,3,3);
    imagesc(xGrid, zGrid, log10(quickDensity + 1));
    axis image; colormap hot; colorbar;
    title(sprintf('Quick density (%d frames, thresh=%d)', numCompound, bestThresh));
    xlabel('Lat [mm]'); ylabel('Ax [mm]');

    sgtitle(sprintf('Detection Tuning (SVD cut=%d, %d frames)', tuneCutoff, numCompound));

catch ME
    fprintf('  SVD section error: %s\n', ME.message);
    fprintf('  %s\n', ME.getReport('basic'));
end

%% Step 6: Noise estimation from compounding strategies
fprintf('\n[6/7] Noise estimation (from compounding comparison)...\n');

if hasPI
    noiseImgs = {img_bmode, img_bmode_blank, img_incoherent, img_piCompound, img_piCompound_blank};
    noiseLabels = {'B-mode (orig)', 'B-mode (blanked)', 'Incoherent (no PI)', ...
                   'PI compound (orig)', 'PI compound (blanked)'};
else
    noiseImgs = {img_single, img_incoherent, img_coh_pos, img_0deg_iq};
    noiseLabels = {'Single 0°', sprintf('Incoherent %d-angle', numAngles), ...
                   sprintf('Coherent %d-angle', numAngles), '0° IQ'};
end
for ni = 1:numel(noiseImgs)
    img = noiseImgs{ni};
    nReg = img(1:20, :);
    snr = 20*log10(max(img(:)) / (std(nReg(:)) + eps));
    fprintf('  %-30s SNR=%.1f dB, max=%.1f, noise std=%.1f\n', ...
        noiseLabels{ni}, snr, max(img(:)), std(nReg(:)));
end

%% Step 7: Data size estimate
fprintf('\n[7/7] Block summary...\n');
nTotalEvt = numel(Config.PulseSequences(1).Events);
bytesPerEvt = size(VadaMode(1).Data,1) * size(VadaMode(1).Data,2) * 2;
nCompFrames = floor(nTotalEvt / eventsPerFrame);
fprintf('  Total events: %d\n', nTotalEvt);
if hasPI
    fprintf('  Events per compound frame: %d (%d angles x 2 PI polarities)\n', eventsPerFrame, numAngles);
else
    fprintf('  Events per compound frame: %d (%d angles x 1 event each, no PI)\n', eventsPerFrame, numAngles);
end
fprintf('  Compound frames: %d\n', nCompFrames);
fprintf('  Transducer: %s (%d elements, pitch=%.4f mm)\n', TxrParam.Name, TxrParam.ArrayNumElements, pitch_mm);
if numel(VadaMode) > eventsPerFrame
    fprintf('  Effective compound frame rate: %.0f Hz\n', ...
        1000 / (VadaMode(eventsPerFrame+1).Timestamp - VadaMode(1).Timestamp));
end
fprintf('  Effective single-event rate: %.0f Hz\n', ...
    1000 / (VadaMode(2).Timestamp - VadaMode(1).Timestamp));
fprintf('  Est. data size: %.1f GB\n', nTotalEvt * bytesPerEvt / 1e9);
fprintf('  Timestamps: %.3f to %.3f ms (loaded)\n', ...
    VadaMode(1).Timestamp, VadaMode(end).Timestamp);

fprintf('\n=== Explore complete ===\n');
if hasPI
    fprintf('If SVD-only images show better bubble contrast than PI images,\n');
    fprintf('use SVD-only mode in the pipeline (skip PI summation).\n');
else
    fprintf('No PI in this sequence. Pipeline will use raw events directly.\n');
    fprintf('SVD clutter filtering will handle tissue suppression.\n');
end
if useGPU, reset(gpuDevice); end
