%% BMODE_RF_EXPLORE.m
% Load F2 B-Mode exported data and check for bubble presence.
%
% Supports two file types:
%   .rf.bmode  -> RF data (has phase, full IQ SVD possible)
%   .raw.bmode -> Envelope-detected data (no phase, power SVD only)
%
% If you only have .raw files, check "Export RF Data" in VevoLab for
% future exports. The .raw data can still confirm bubble presence via
% temporal variance and power SVD.
%
% Author: Eli Wirth-Apley / Sun Lab, Northeastern University, February 2026

clearvars; close all; clc;

%% CONFIGURATION
dataFolder   = 'C:\path\to\VADA_data';  % CHANGE
baseFilename = '';  % CHANGE (no extension)
modeName     = '.bmode';
vadaScriptsPath = 'C:\path\to\VevoF2\RFExampleScriptsF2\Matlab';
useGPU       = true;

% Data type: 'rf' or 'raw' (auto-detected if left empty)
dataType     = '';  % Leave empty to auto-detect

addpath(genpath(vadaScriptsPath));

%% Step 1: Detect file type and load
fprintf('=== B-Mode Explorer ===\n');

% Auto-detect file type
if isempty(dataType)
    rfExists  = exist(fullfile(dataFolder, [baseFilename '.rf' modeName]), 'file');
    rawExists = exist(fullfile(dataFolder, [baseFilename '.raw' modeName]), 'file');
    if rfExists
        dataType = 'rf';
    elseif rawExists
        dataType = 'raw';
    else
        error('No .rf%s or .raw%s file found for %s', modeName, modeName, baseFilename);
    end
end
fprintf('  Data type: %s\n', dataType);

switch dataType
    case 'rf'
        fprintf('  RF data: full complex IQ SVD available\n');
        hasPhase = true;
        
        fprintf('[1/5] Loading RF data...\n');
        BMode = VsiBModeRfRead(dataFolder, baseFilename, -1, modeName);
        numFrames = numel(BMode.Data);
        nZ = numel(BMode.Depth);
        nX = numel(BMode.Width);
        
        % Timestamps: RF has per-line timestamps
        if size(BMode.Timestamp, 2) > 1
            ts_ms = BMode.Timestamp(:, 1) * 1000;
        else
            ts_ms = BMode.Timestamp(:) * 1000;
        end
        
    case 'raw'
        fprintf('  RAW data: envelope only (no phase). Power SVD + temporal variance.\n');
        fprintf('  Tip: export with "RF Data" checked for full IQ SVD next time.\n');
        hasPhase = false;
        
        fprintf('[1/5] Loading RAW data...\n');
        BMode = VsiBModeRawRead(dataFolder, baseFilename, -1, modeName);
        numFrames = numel(BMode.Data);
        nZ = numel(BMode.Depth);
        nX = numel(BMode.Width);
        ts_ms = BMode.Timestamp(:);
        
    otherwise
        error('dataType must be ''rf'' or ''raw''');
end

fprintf('  Array: %s, %d elements, pitch=%.3f mm\n', ...
    BMode.ArrayType, BMode.NumElem, BMode.ElemPitch);
fprintf('  Depth: %.1f to %.1f mm (%d samples)\n', ...
    min(BMode.Depth), max(BMode.Depth), nZ);
fprintf('  Width: %.1f to %.1f mm (%d lines)\n', ...
    min(BMode.Width), max(BMode.Width), nX);
fprintf('  Frames: %d\n', numFrames);

if numFrames > 1
    dt_ms = median(diff(ts_ms));
    frameRate = 1000 / dt_ms;
    fprintf('  Frame rate: %.1f Hz (dt=%.2f ms)\n', frameRate, dt_ms);
else
    frameRate = 20;
    fprintf('  Single frame only.\n');
end

%% Step 2: Build image stack
fprintf('\n[2/5] Building image stack...\n');

if hasPhase
    % RF data: Hilbert transform for complex IQ
    iqStack = zeros(nZ, nX, numFrames, 'single');
    envStack = zeros(nZ, nX, numFrames, 'single');
    for f = 1:numFrames
        iq = single(hilbert(double(BMode.Data{f})));
        iqStack(:,:,f) = iq;
        envStack(:,:,f) = abs(iq);
    end
else
    % RAW data: already envelope-detected
    envStack = zeros(nZ, nX, numFrames, 'single');
    for f = 1:numFrames
        envStack(:,:,f) = single(BMode.Data{f});
    end
end

%% Step 3: Display B-Mode + temporal analysis
fprintf('\n[3/5] B-Mode display and temporal analysis...\n');

figure('Name', 'B-Mode Overview', 'Position', [100 100 1400 500]);

% Single frame
subplot(1,4,1);
img1 = envStack(:,:,1);
imagesc(BMode.Width, BMode.Depth, 20*log10(img1/max(img1(:))+eps));
axis image; colormap gray; colorbar; caxis([-60 0]);
xlabel('Width [mm]'); ylabel('Depth [mm]');
title('Frame 1');

% Mean
subplot(1,4,2);
meanEnv = mean(envStack, 3);
imagesc(BMode.Width, BMode.Depth, 20*log10(meanEnv/max(meanEnv(:))+eps));
axis image; colormap gray; colorbar; caxis([-60 0]);
xlabel('Width [mm]'); ylabel('Depth [mm]');
title(sprintf('Mean (%d fr)', numFrames));

% Temporal std (KEY: bubbles = high variance)
subplot(1,4,3);
stdImg = std(envStack, 0, 3);
imagesc(BMode.Width, BMode.Depth, 20*log10(stdImg/max(stdImg(:))+eps));
axis image; colormap gray; colorbar; caxis([-60 0]);
xlabel('Width [mm]'); ylabel('Depth [mm]');
title('Temporal std (MOTION)');

% Temporal std normalized by mean (coefficient of variation)
% This highlights regions where variance is high relative to signal level
subplot(1,4,4);
cv = stdImg ./ (meanEnv + eps);
imagesc(BMode.Width, BMode.Depth, cv);
axis image; colormap hot; colorbar;
caxis([0, prctile(cv(:), 99)]);
xlabel('Width [mm]'); ylabel('Depth [mm]');
title('Coeff of Variation (std/mean)');

sgtitle(sprintf('B-Mode %s: %d frames at %.0f Hz', upper(dataType), numFrames, frameRate));

%% Step 4: SVD (complex IQ if RF, power if RAW)
fprintf('\n[4/5] SVD analysis...\n');

if hasPhase
    fprintf('  Using complex IQ SVD (RF data)\n');
    svdInput = reshape(iqStack, nZ*nX, numFrames);
else
    fprintf('  Using power SVD (envelope data, no phase)\n');
    % For envelope data, SVD can still separate stationary (tissue)
    % from fluctuating (bubble) components based on temporal intensity changes
    svdInput = reshape(envStack, nZ*nX, numFrames);
end

if useGPU && gpuDeviceCount > 0
    [U, S, V] = svd(gpuArray(single(svdInput)), 'econ');
    U = gather(U); S = gather(S); V = gather(V);
else
    [U, S, V] = svd(single(svdInput), 'econ');
end
singVals = diag(S);
clear svdInput;

fprintf('  SV1/SV2: %.1f, SV1/SV5: %.1f\n', ...
    singVals(1)/singVals(2), singVals(1)/singVals(min(5,end)));

% Display SVD filtered at multiple cutoffs
svdCutoffs = [1, 2, 3, 5, 10];
svdCutoffs = svdCutoffs(svdCutoffs < numFrames);

figure('Name', 'SVD Clutter Filter', 'Position', [50 50 1800 700]);
nCols = ceil((numel(svdCutoffs)+2) / 2);

subplot(2, nCols, 1);
imagesc(BMode.Width, BMode.Depth, 20*log10(meanEnv/max(meanEnv(:))+eps));
axis image; colormap gray; colorbar; caxis([-60 0]);
title('Mean (tissue)'); xlabel('W [mm]'); ylabel('D [mm]');

subplot(2, nCols, 2);
if hasPhase
    cStd = abs(std(iqStack, 0, 3));
else
    cStd = std(envStack, 0, 3);
end
imagesc(BMode.Width, BMode.Depth, 20*log10(cStd/max(cStd(:))+eps));
axis image; colormap gray; colorbar; caxis([-60 0]);
title('Temporal std'); xlabel('W [mm]'); ylabel('D [mm]');

for ci = 1:numel(svdCutoffs)
    cutoff = svdCutoffs(ci);
    nKeep = size(U, 2);
    if cutoff < nKeep
        filtered = U(:, cutoff+1:nKeep) * S(cutoff+1:nKeep, cutoff+1:nKeep) * V(:, cutoff+1:nKeep)';
    else
        filtered = zeros(nZ*nX, numFrames, 'single');
    end
    filtStack = reshape(filtered, nZ, nX, numFrames);
    filtMean = mean(abs(filtStack), 3);
    
    subplot(2, nCols, ci+2);
    imagesc(BMode.Width, BMode.Depth, 20*log10(filtMean/max(filtMean(:))+eps));
    axis image; colormap gray; colorbar; caxis([-60 0]);
    title(sprintf('SVD cut=%d', cutoff));
    xlabel('W [mm]'); ylabel('D [mm]');
end
sgtitle(sprintf('SVD Filter (%s, %d frames)', upper(dataType), numFrames));

% Individual frames
figure('Name', 'Individual Filtered Frames', 'Position', [50 50 1600 800]);
for tc = 1:min(2, numel(svdCutoffs))
    testCut = svdCutoffs(min(tc*2, numel(svdCutoffs)));
    if testCut >= size(U, 2), continue; end
    filt = U(:, testCut+1:end) * S(testCut+1:end, testCut+1:end) * V(:, testCut+1:end)';
    fFrames = reshape(filt, nZ, nX, numFrames);
    
    frameIdxs = round(linspace(1, numFrames, 4));
    rowOff = (tc-1) * 4;
    for fi = 1:4
        subplot(2, 4, fi + rowOff);
        frame = abs(fFrames(:,:,frameIdxs(fi)));
        imagesc(BMode.Width, BMode.Depth, frame);
        axis image; colormap hot; colorbar;
        title(sprintf('Fr %d (cut=%d)', frameIdxs(fi), testCut), 'FontSize', 9);
        xlabel('W [mm]'); ylabel('D [mm]');
    end
end
sgtitle('Individual filtered frames');

%% Step 5: SV spectrum
fprintf('\n[5/5] SV spectrum...\n');

figure('Name', 'SV Spectrum', 'Position', [100 100 700 400]);
nShow = min(numFrames, 30);
subplot(1,2,1);
semilogy(singVals/singVals(1), 'b.-', 'LineWidth', 1.5);
xlabel('SV Index'); ylabel('Normalized'); title('Full spectrum'); grid on;
hold on;
for ci = 1:numel(svdCutoffs)
    xline(svdCutoffs(ci)+0.5, '--', sprintf('%d', svdCutoffs(ci)));
end

subplot(1,2,2);
bar(singVals(1:nShow)/singVals(1), 'FaceColor', [0.3 0.5 0.8]);
xlabel('SV Index'); ylabel('Normalized'); title(sprintf('First %d SVs', nShow)); grid on;
sgtitle(sprintf('SVD Spectrum (%s, %d frames)', upper(dataType), numFrames));

% SNR
fprintf('\n  SNR by cutoff:\n');
for testCut = svdCutoffs
    if testCut >= size(U, 2), continue; end
    filt = U(:, testCut+1:end) * S(testCut+1:end, testCut+1:end) * V(:, testCut+1:end)';
    fImg = mean(abs(reshape(filt, nZ, nX, numFrames)), 3);
    nReg = fImg(1:min(20,nZ), :);
    snr = 20*log10(max(fImg(:)) / (std(nReg(:)) + eps));
    fprintf('    cut=%2d: SNR=%.1f dB\n', testCut, snr);
end

fprintf('\n=== B-Mode Explore complete ===\n');
if ~hasPhase
    fprintf('\nNOTE: This was RAW (envelope) data. SVD is power-based only.\n');
    fprintf('For full complex IQ SVD, re-export from VevoLab with "RF Data" checked.\n');
end
