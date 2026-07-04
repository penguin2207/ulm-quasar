%% PI_DEBUG_COMPOUND.m
% Comprehensive PI quality debugging and coherent compound imaging test.
%
% Tests PI cancellation with and without RF blanking for steered angles,
% and compares beamformed images across compounding strategies.
%
% The key fix: at steered angles, the first N RF samples are contaminated
% by TX cross-talk during the element-sequential firing window. N is:
%   N = (nRx-1) * pitch * |sin(theta)| / c * fs
% Blanking these samples before PI summation improves cancellation.
%
% Author: Eli Wirth-Apley / Sun Lab, Northeastern University, February 2026

clearvars; close all; clc;

%% CONFIGURATION
dataFolder      = 'C:\path\to\VADA_data';
baseFilename    = '';  % CHANGE (no extension)
modeName        = '.vada';
vadaScriptsPath = 'C:\path\to\VevoF2\RFExampleScriptsF2\Matlab';
sosOverride     = [];    % Set to 1540 for agarose, [] for XML value
useGPU          = true;
blankMargin     = 1.5;  % Base safety margin (auto-scaled by voltage)
voltageOverride = [];   % TX voltage override [%]. [] = read from XML.

addpath(genpath(vadaScriptsPath));
if useGPU && gpuDeviceCount > 0
    g = gpuDevice; reset(g);
    fprintf('GPU: %s (%.1f GB VRAM)\n', g.Name, g.TotalMemory/1e9);
end

%% STEP 1: Load first 3 compound frames
fprintf('\n[1/6] Loading data...\n');
numLoad = 18;  % 3 compound frames x 6 events
[VadaMode, Param, TxrParam, Config] = VsiVadaDataRead(dataFolder, baseFilename, 1:numLoad, modeName);

% Pitch detection
rawPitch = TxrParam.ArrayPitch;
if rawPitch == 0 || isnan(rawPitch), pitch_mm = 0.300;
elseif rawPitch < 10, pitch_mm = rawPitch;
else, pitch_mm = rawPitch / 1000; end

c = Param.SoSMedia; if c == 0, c = 1540; end
if ~isempty(sosOverride), c = sosOverride; end
fs_MHz = Param.SampleFreq;
depthOffset_mm = Param.DepthOffset;
nRx = numel(VadaMode(1).Elements);

% Read TX voltage from XML
if isempty(voltageOverride)
    xmlPath = fullfile(dataFolder, [baseFilename modeName '.xml']);
    if exist(xmlPath, 'file')
        xmlText = fileread(xmlPath);
        tokens = regexp(xmlText, '<parameter\s+name="([^"]+)"\s+value="([^"]*)"', 'tokens');
        xmlP = struct();
        for ti = 1:numel(tokens)
            fn = strrep(strrep(tokens{ti}{1}, '-', '_'), '/', '_');
            xmlP.(fn) = tokens{ti}{2};
        end
        if isfield(xmlP, 'Vada_Mode_Voltage_Rail_High')
            txVoltage = str2double(xmlP.Vada_Mode_Voltage_Rail_High);
            fprintf('  TX voltage: %s-%s%% (from XML)\n', ...
                xmlP.Vada_Mode_Voltage_Rail_Low, xmlP.Vada_Mode_Voltage_Rail_High);
        else
            txVoltage = [];
            fprintf('  TX voltage: not found in XML\n');
        end
    else
        txVoltage = [];
    end
else
    txVoltage = voltageOverride;
    fprintf('  TX voltage: %.0f%% (manual override)\n', txVoltage);
end

elemPos_mm = ((1:TxrParam.ArrayNumElements) - (TxrParam.ArrayNumElements+1)/2) * pitch_mm;
rxElem = VadaMode(1).Elements;
rxPos_mm = elemPos_mm(rxElem);

fprintf('  %s | %d RX elem | pitch=%.3f mm | Fs=%.1f MHz | SoS=%.0f m/s\n', ...
    TxrParam.Name, nRx, pitch_mm, fs_MHz, c);

%% STEP 2: Waveform analysis
fprintf('\n[2/6] Waveform analysis...\n');

for ev = 1:6
    wf = VadaMode(ev).Waveform;
    if isempty(wf.Channel), fprintf('  Event %d: NO WAVEFORM\n', ev); continue; end
    ch = wf.Channel(1);
    fprintf('  Event %d: angle=%+.1f, invert=%d, freq=%.1fMHz, samples=[%s]\n', ...
        ev, VadaMode(ev).TxDelay.angle, ch.invert, ch.frequency, num2str(ch.samples(:)'));
end

fprintf('\n  Waveform pair comparison:\n');
for pair = 1:3
    posIdx = (pair-1)*2 + 1;
    negIdx = (pair-1)*2 + 2;
    wfPos = VadaMode(posIdx).Waveform.Channel(1).samples;
    wfNeg = VadaMode(negIdx).Waveform.Channel(1).samples;
    
    fprintf('  Pair %d (%+.1f°): exact negation=%s, sum=[%s]\n', ...
        pair, VadaMode(posIdx).TxDelay.angle, ...
        string(all(wfPos == -wfNeg)), num2str((wfPos+wfNeg)'));
end

%% STEP 3: Compute blanking parameters
fprintf('\n[3/6] Steering blanking parameters...\n');

angles = zeros(3,1);
blankInfo = struct();
for pair = 1:3
    angles(pair) = VadaMode((pair-1)*2+1).TxDelay.angle;
    bi = compute_steering_blanking(angles(pair), nRx, pitch_mm, c, fs_MHz, blankMargin, txVoltage);
    blankInfo(pair).nBlank = bi.nBlank;
    blankInfo(pair).minDepth_mm = bi.minDepth_mm;
    blankInfo(pair).delaySpread_us = bi.delaySpread_us;
    blankInfo(pair).delaySpread_samples = bi.delaySpread_samples;
    
    fprintf('  Angle %+5.1f: delay spread=%.2f us (%.0f samples), blank=%d samples, minDepth=%.2f mm\n', ...
        angles(pair), bi.delaySpread_us, bi.delaySpread_samples, bi.nBlank, bi.minDepth_mm);
end

%% STEP 4: PI cancellation comparison (with and without blanking)
fprintf('\n[4/6] PI cancellation: original vs blanked...\n');

figure('Name', 'PI Cancellation: Original vs Blanked', 'Position', [50 50 1800 900]);

for pair = 1:3
    posIdx = (pair-1)*2 + 1;
    negIdx = (pair-1)*2 + 2;
    
    rfPos = double(VadaMode(posIdx).Data);
    rfNeg = double(VadaMode(negIdx).Data);
    nSamp = size(rfPos, 1);
    nCh = size(rfPos, 2);
    nb = blankInfo(pair).nBlank;
    
    % --- Original (no blanking) ---
    rfSum_orig = rfPos + rfNeg;
    
    % Per-channel cancellation (original)
    cancOrig = zeros(nCh, 1);
    for ch = 1:nCh
        cancOrig(ch) = 10*log10(sum(rfSum_orig(:,ch).^2) / (sum(rfPos(:,ch).^2)+eps) + eps);
    end
    
    % --- Blanked ---
    rfPos_b = rfPos; rfNeg_b = rfNeg;
    if nb > 0 && nb < nSamp
        rfPos_b(1:nb, :) = 0;
        rfNeg_b(1:nb, :) = 0;
    end
    rfSum_blank = rfPos_b + rfNeg_b;
    
    % Per-channel cancellation (blanked)
    cancBlank = zeros(nCh, 1);
    for ch = 1:nCh
        cancBlank(ch) = 10*log10(sum(rfSum_blank(:,ch).^2) / (sum(rfPos_b(:,ch).^2)+eps) + eps);
    end
    
    fprintf('  Pair %d (%+.1f°): blank=%d samples\n', pair, angles(pair), nb);
    fprintf('    Original: mean=%.1f dB, range=[%.1f to %.1f]\n', ...
        mean(cancOrig), min(cancOrig), max(cancOrig));
    fprintf('    Blanked:  mean=%.1f dB, range=[%.1f to %.1f]\n', ...
        mean(cancBlank), min(cancBlank), max(cancBlank));
    fprintf('    Improvement: %.1f dB\n', mean(cancOrig) - mean(cancBlank));
    
    % Plot: per-channel cancellation comparison
    subplot(3, 4, (pair-1)*4 + 1);
    plot(cancOrig, 'r.-', 'DisplayName', 'Original'); hold on;
    plot(cancBlank, 'b.-', 'DisplayName', sprintf('Blanked (%d samp)', nb));
    legend('Location', 'best'); xlabel('RX Channel'); ylabel('Cancel [dB]');
    title(sprintf('%+.1f°: per-channel', angles(pair)));
    yline(-10, 'g--'); yline(-20, 'g--');
    
    % Plot: RF overlay (mid channel)
    midCh = round(nCh/2);
    subplot(3, 4, (pair-1)*4 + 2);
    plot(rfPos(:,midCh), 'b'); hold on;
    plot(-rfNeg(:,midCh), 'r--');
    if nb > 0, xline(nb, 'k-', sprintf('blank=%d', nb), 'LineWidth', 2); end
    title(sprintf('%+.1f°: RF +/-(neg) ch%d', angles(pair), midCh));
    xlabel('Sample'); xlim([1 min(200, nSamp)]);
    
    % Plot: windowed cancellation map (original)
    winLen = 50;
    nWin = floor(nSamp/winLen);
    cancMapOrig = zeros(nWin, nCh);
    cancMapBlank = zeros(nWin, nCh);
    for w = 1:nWin
        idx = (w-1)*winLen+1 : w*winLen;
        cancMapOrig(w,:) = 10*log10(sum(rfSum_orig(idx,:).^2,1) ./ (sum(rfPos(idx,:).^2,1)+eps)+eps);
        cancMapBlank(w,:) = 10*log10(sum(rfSum_blank(idx,:).^2,1) ./ (sum(rfPos_b(idx,:).^2,1)+eps)+eps);
    end
    
    subplot(3, 4, (pair-1)*4 + 3);
    imagesc(1:nCh, (1:nWin)*winLen, cancMapOrig);
    colorbar; caxis([-20 5]); colormap(gca, 'jet');
    xlabel('Channel'); ylabel('Sample');
    title(sprintf('%+.1f° cancel map (ORIGINAL)', angles(pair)));
    if nb > 0, yline(nb, 'k-', 'LineWidth', 2); end
    
    subplot(3, 4, (pair-1)*4 + 4);
    imagesc(1:nCh, (1:nWin)*winLen, cancMapBlank);
    colorbar; caxis([-20 5]); colormap(gca, 'jet');
    xlabel('Channel'); ylabel('Sample');
    title(sprintf('%+.1f° cancel map (BLANKED)', angles(pair)));
    if nb > 0, yline(nb, 'k-', 'LineWidth', 2); end
end
sgtitle('PI Cancellation: Original vs RF Blanking');

%% STEP 5: Beamformed compound comparison (with and without blanking)
fprintf('\n[5/6] Beamforming comparison...\n');

xGrid = -10:0.050:10;
zGrid = 2:0.025:25;

% Precompute delay tables
dtables = cell(3,1);
for a = 1:3
    dtables{a} = beamform_planewave_gpu([], rxPos_mm, angles(a), ...
        xGrid, zGrid, fs_MHz, c, depthOffset_mm, []);
end

% Beamform all 6 events
bfEvents = zeros(numel(zGrid), numel(xGrid), 6);
wb = waitbar(0, 'Beamforming...', 'Name', 'PI Debug');
for ev = 1:6
    waitbar(ev/6, wb);
    angleIdx = ceil(ev/2);
    bfEvents(:,:,ev) = beamform_planewave_gpu(single(VadaMode(ev).Data), rxPos_mm, ...
        angles(angleIdx), xGrid, zGrid, fs_MHz, c, depthOffset_mm, dtables{angleIdx});
end
close(wb);

% Beamform blanked PI at each angle
bfPI_orig = zeros(numel(zGrid), numel(xGrid), 3);
bfPI_blank = zeros(numel(zGrid), numel(xGrid), 3);
bfBmode_orig = zeros(numel(zGrid), numel(xGrid), 3);
bfBmode_blank = zeros(numel(zGrid), numel(xGrid), 3);

for pair = 1:3
    posIdx = (pair-1)*2 + 1;
    negIdx = (pair-1)*2 + 2;
    nb = blankInfo(pair).nBlank;
    
    % Original PI
    rfPI_orig = single(VadaMode(posIdx).Data) + single(VadaMode(negIdx).Data);
    rfBm_orig = single(VadaMode(posIdx).Data) - single(VadaMode(negIdx).Data);
    
    % Blanked PI
    rfPos = single(VadaMode(posIdx).Data);
    rfNeg = single(VadaMode(negIdx).Data);
    if nb > 0
        rfPos(1:nb,:) = 0;
        rfNeg(1:nb,:) = 0;
    end
    rfPI_blank = rfPos + rfNeg;
    rfBm_blank = rfPos - rfNeg;
    
    bfPI_orig(:,:,pair) = beamform_planewave_gpu(rfPI_orig, rxPos_mm, angles(pair), ...
        xGrid, zGrid, fs_MHz, c, depthOffset_mm, dtables{pair});
    bfPI_blank(:,:,pair) = beamform_planewave_gpu(rfPI_blank, rxPos_mm, angles(pair), ...
        xGrid, zGrid, fs_MHz, c, depthOffset_mm, dtables{pair});
    bfBmode_orig(:,:,pair) = beamform_planewave_gpu(rfBm_orig, rxPos_mm, angles(pair), ...
        xGrid, zGrid, fs_MHz, c, depthOffset_mm, dtables{pair});
    bfBmode_blank(:,:,pair) = beamform_planewave_gpu(rfBm_blank, rxPos_mm, angles(pair), ...
        xGrid, zGrid, fs_MHz, c, depthOffset_mm, dtables{pair});
end

% Compound images
img_0deg_only = abs(bfEvents(:,:,3));  % 0deg positive event
img_incoh = abs(bfEvents(:,:,1)) + abs(bfEvents(:,:,3)) + abs(bfEvents(:,:,5));

img_compound_pi_orig = abs(sum(bfPI_orig, 3));
img_compound_pi_blank = abs(sum(bfPI_blank, 3));
img_compound_bm_orig = abs(sum(bfBmode_orig, 3));
img_compound_bm_blank = abs(sum(bfBmode_blank, 3));

% Also try coherent compound without PI (all 6 events summed)
img_coh_all = abs(sum(bfEvents, 3));
% Coherent compound positive events only (3 events, no PI)
img_coh_pos = abs(bfEvents(:,:,1) + bfEvents(:,:,3) + bfEvents(:,:,5));

figure('Name', 'Beamforming Comparison', 'Position', [50 50 1800 800]);
titles = {'0° only', 'Incoherent 3-angle', 'Coherent 3-pos (no PI)', ...
          'PI compound (ORIGINAL)', 'PI compound (BLANKED)', ...
          'B-mode compound (ORIG)', 'B-mode compound (BLANK)', ...
          'Coherent all 6'};
imgs = {img_0deg_only, img_incoh, img_coh_pos, ...
        img_compound_pi_orig, img_compound_pi_blank, ...
        img_compound_bm_orig, img_compound_bm_blank, ...
        img_coh_all};

for i = 1:8
    subplot(2, 4, i);
    img = imgs{i};
    imagesc(xGrid, zGrid, 20*log10(img / max(img(:)) + eps));
    axis image; colormap gray; colorbar; caxis([-60 0]);
    xlabel('Lat [mm]'); ylabel('Ax [mm]'); title(titles{i}, 'FontSize', 9);
end
sgtitle('Compound Imaging: Effect of RF Blanking on PI Quality');

%% STEP 6: Cross-correlation lag analysis
fprintf('\n[6/6] Cross-correlation lag analysis...\n');

figure('Name', 'PI Time Alignment', 'Position', [100 100 1200 400]);
colors = {'b', 'r', 'g'};
hold on;

for pair = 1:3
    posIdx = (pair-1)*2 + 1;
    negIdx = (pair-1)*2 + 2;
    nCh = size(VadaMode(posIdx).Data, 2);
    lagPerCh = zeros(nCh, 1);
    for ch = 1:nCh
        [xc, lags] = xcorr(double(VadaMode(posIdx).Data(:,ch)), ...
                           double(VadaMode(negIdx).Data(:,ch)), 5);
        [~, pk] = max(abs(xc));
        lagPerCh(ch) = lags(pk);
    end
    plot(lagPerCh, '.-', 'Color', colors{pair}, 'DisplayName', ...
        sprintf('%+.1f° (mean lag=%.2f)', angles(pair), mean(lagPerCh)));
end
legend('Location', 'best'); xlabel('RX Channel'); ylabel('Lag (samples)');
title('PI Time Alignment: lag per channel per angle');
yline(0, 'k--'); grid on;

%% Summary
fprintf('\n=== PI DEBUG SUMMARY ===\n');
fprintf('Blanking margin: %.1fx\n', blankMargin);
fprintf('\n  %-8s %8s %12s %12s %10s\n', 'Angle', 'nBlank', 'Orig [dB]', 'Blank [dB]', 'Improve');
fprintf('  %s\n', repmat('-', 1, 52));
for pair = 1:3
    posIdx = (pair-1)*2 + 1;
    rfPos = double(VadaMode(posIdx).Data);
    rfNeg = double(VadaMode((pair-1)*2+2).Data);
    nb = blankInfo(pair).nBlank;
    
    % Original
    s = rfPos + rfNeg;
    co = 10*log10(mean(sum(s.^2,1)) / (mean(sum(rfPos.^2,1))+eps) + eps);
    
    % Blanked
    rp = rfPos; rn = rfNeg;
    if nb > 0, rp(1:nb,:)=0; rn(1:nb,:)=0; end
    sb = rp + rn;
    cb = 10*log10(mean(sum(sb.^2,1)) / (mean(sum(rp.^2,1))+eps) + eps);
    
    fprintf('  %+5.1f°   %8d %11.1f %11.1f %9.1f\n', ...
        angles(pair), nb, co, cb, co - cb);
end

fprintf('\nRecommendation:\n');
fprintf('  If blanked cancellation < -10 dB at all angles:\n');
fprintf('    -> Use config.zeroOnly=false with config.blankSteering=true\n');
fprintf('  If blanked cancellation still poor at steered angles:\n');
fprintf('    -> Use config.zeroOnly=true (0deg only, SVD handles tissue)\n');
fprintf('\n=== Done ===\n');
