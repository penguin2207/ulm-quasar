%% INSPECT_SVD_FRAMES.m
% Memory-efficient SVD-filtered frame viewer with video export.
%
% Loads cached IQ data from a PARAMETER_TUNER or BATCH_PARAMETER_SWEEP
% cache file and visualizes SVD-filtered frames to distinguish real
% microbubble signal from tissue residual.
%
% Uses truncated SVD (svds) to compute ONLY the components being removed,
% then subtracts them frame-by-frame — much more memory-efficient than
% full SVD for large frame counts.
%
% OUTPUTS (saved to outputFolder):
%   1. Single-panel video: SVD-filtered frames (MP4)
%   2. Three-panel video:  current frame | temporal MIP | running accumulator
%   3. Accumulator snapshot: final max-projection across all frames (PNG)
%
% The temporal MIP (panel 2) shows the max over a sliding window — real
% bubbles show as short streaks that move; tissue residual shows as
% persistent blobs in the same spot.
%
% Author: Eli Wirth-Apley / Sun Lab, Northeastern University, March 2026

clearvars; close all; clc;

%% === CONFIGURATION ===

% --- Input: path to cached IQ .mat file ---
% From PARAMETER_TUNER: contains IQ_raw, xGrid, zGrid
% From BATCH_PARAMETER_SWEEP: block_XX_IQ.mat contains IQ_raw
cacheFile = '';  % <-- SET THIS to your cached IQ .mat file path

% --- SVD settings ---
svdCut      = 2;       % SVD components to remove (try 1, 5, 10, 20, 50)
maxFrames   = 0;       % 0 = use all frames, otherwise subsample to this count
useGPU      = true;    % GPU for SVD
contrastPct = 99.5;    % Percentile for display clipping

% --- Video settings ---
saveVideo    = true;    % true = save MP4s + PNG, false = live playback only
videoFPS     = 15;      % Playback framerate for saved video
playSpeed    = 0.05;    % Seconds per frame in live mode (0 = click-to-advance)

% --- Temporal MIP window ---
mipWindow = 5;          % Frames to max-project over (shows short bubble streaks)

% --- Output folder ('' = same folder as cacheFile) ---
outputFolder = '';

%% === SETUP ===
if isempty(cacheFile)
    [f, p] = uigetfile('*.mat', 'Select cached IQ .mat file');
    if isequal(f, 0), error('No file selected.'); end
    cacheFile = fullfile(p, f);
end

fprintf('Loading: %s\n', cacheFile);
d = load(cacheFile);

% Handle different cache formats
if isfield(d, 'IQ_raw')
    IQ = d.IQ_raw;
elseif isfield(d, 'IQ_block')
    IQ = d.IQ_block;
else
    fn = fieldnames(d);
    % Find first 3D complex array
    IQ = [];
    for i = 1:numel(fn)
        v = d.(fn{i});
        if ndims(v) == 3 && ~isreal(v)
            IQ = v;
            fprintf('  Using variable: %s\n', fn{i});
            break;
        end
    end
    if isempty(IQ), error('No 3D complex IQ array found in cache file.'); end
end

% Grid axes
if isfield(d, 'xGrid')
    xGrid = d.xGrid;
    zGrid = d.zGrid;
else
    % Fallback: pixel indices
    [nZ, nX, ~] = size(IQ);
    xGrid = 1:nX;
    zGrid = 1:nZ;
    fprintf('  WARNING: No grid axes found, using pixel indices.\n');
end

[nZ, nX, nF] = size(IQ);

% Subsample frames if requested
if maxFrames > 0 && nF > maxFrames
    frameIdx = round(linspace(1, nF, maxFrames));
    IQ = IQ(:,:,frameIdx);
    nF = maxFrames;
    fprintf('Subsampled to %d frames\n', nF);
end

nPix = nZ * nX;
fprintf('Data: %d x %d x %d frames (%.1f MB)\n', nZ, nX, nF, nPix*nF*8/1e6);

if isempty(outputFolder)
    outputFolder = fileparts(cacheFile);
end
if saveVideo && ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

%% === TRUNCATED SVD (memory-efficient) ===
fprintf('Computing truncated SVD (first %d components)...\n', svdCut);
tSVD = tic;

M = reshape(IQ, nPix, nF);
clear IQ;

if useGPU && gpuDeviceCount > 0
    M = gpuArray(M);
end

[U_k, S_k, V_k] = svds(double(M), svdCut);
US_k = U_k * S_k;
clear U_k S_k;

if useGPU && gpuDeviceCount > 0
    US_k = gather(US_k);
    V_k  = gather(V_k);
    M    = gather(M);
end

fprintf('SVD done in %.1f s. Tissue subspace: %d components.\n', toc(tSVD), svdCut);

%% === COMPUTE GLOBAL CONTRAST ===
fprintf('Computing display range...\n');
sampleIdx = round(linspace(1, nF, min(30, nF)));
sampVals = zeros(numel(sampleIdx), 1);
for i = 1:numel(sampleIdx)
    f = sampleIdx(i);
    fr = abs(single(M(:,f) - US_k * V_k(f,:)'));
    sampVals(i) = prctile(fr, contrastPct);
end
clim_val = median(sampVals);
fprintf('Display clim: [0, %.2f]\n', clim_val);

%% === VIDEO 1: Single-panel SVD-filtered frames ===
if saveVideo
    vidFile1 = fullfile(outputFolder, sprintf('svd_cut%d_frames.mp4', svdCut));
    fprintf('\nWriting single-panel video: %s\n', vidFile1);
    vw1 = VideoWriter(vidFile1, 'MPEG-4');
    vw1.FrameRate = videoFPS; vw1.Quality = 95; open(vw1);
end

fig1 = figure('Position', [100 100 960 720], 'Color', 'k');
ax1 = axes(fig1);

fprintf('\nPlaying %d frames (SVD cutoff = %d)...\n', nF, svdCut);
if ~saveVideo && playSpeed == 0
    fprintf('  Click figure to advance. Close figure to stop.\n');
end

for f = 1:nF
    if ~isvalid(fig1), break; end

    frame_2d = reshape(abs(single(M(:,f) - US_k * V_k(f,:)')), nZ, nX);

    imagesc(ax1, xGrid, zGrid, frame_2d);
    colormap(ax1, 'hot'); caxis(ax1, [0 clim_val]); colorbar(ax1);
    xlabel(ax1, 'Lateral [mm]'); ylabel(ax1, 'Axial [mm]');
    title(ax1, sprintf('Frame %d / %d   |   SVD cutoff = %d', f, nF, svdCut), ...
        'Color', 'w', 'FontSize', 14);
    axis(ax1, 'image');
    drawnow;

    if saveVideo
        writeVideo(vw1, getframe(fig1));
        if mod(f, 100) == 0
            fprintf('  %d / %d (%.0f%%)\n', f, nF, f/nF*100);
        end
    else
        if playSpeed == 0
            try waitforbuttonpress; catch, break; end
        else
            pause(playSpeed);
        end
    end
end

if saveVideo
    close(vw1); close(fig1);
    fprintf('  Saved: %s\n', vidFile1);
else
    fprintf('Done. Figure left open for inspection.\n');
end

%% === VIDEO 2: Three-panel + accumulator snapshot ===
if saveVideo
    vidFile2 = fullfile(outputFolder, sprintf('svd_cut%d_3panel.mp4', svdCut));
    fprintf('\nWriting 3-panel video: %s\n', vidFile2);

    fig2 = figure('Position', [50 50 1800 600], 'Color', 'k', 'Visible', 'off');
    vw2 = VideoWriter(vidFile2, 'MPEG-4');
    vw2.FrameRate = videoFPS; vw2.Quality = 95; open(vw2);

    mipBuf = zeros(nZ, nX, mipWindow, 'single');
    accumMax = zeros(nZ, nX, 'single');

    for f = 1:nF
        frame_2d = reshape(abs(single(M(:,f) - US_k * V_k(f,:)')), nZ, nX);

        bufIdx = mod(f-1, mipWindow) + 1;
        mipBuf(:,:,bufIdx) = frame_2d;
        if f >= mipWindow
            tempMIP = max(mipBuf, [], 3);
        else
            tempMIP = max(mipBuf(:,:,1:min(f, mipWindow)), [], 3);
        end
        accumMax = max(accumMax, frame_2d);

        ax2a = subplot(1,3,1, 'Parent', fig2);
        imagesc(ax2a, xGrid, zGrid, frame_2d);
        colormap(ax2a, 'hot'); caxis(ax2a, [0 clim_val]);
        title(ax2a, sprintf('Current frame %d/%d', f, nF), 'Color', 'w');
        ylabel(ax2a, 'Axial [mm]'); xlabel(ax2a, 'Lat [mm]');
        axis(ax2a, 'image');

        ax2b = subplot(1,3,2, 'Parent', fig2);
        imagesc(ax2b, xGrid, zGrid, tempMIP);
        colormap(ax2b, 'hot'); caxis(ax2b, [0 clim_val * 1.5]);
        title(ax2b, sprintf('Temporal MIP (%d-frame window)', mipWindow), 'Color', 'w');
        xlabel(ax2b, 'Lat [mm]');
        axis(ax2b, 'image');

        ax2c = subplot(1,3,3, 'Parent', fig2);
        imagesc(ax2c, xGrid, zGrid, log10(accumMax + 1));
        colormap(ax2c, 'hot'); caxis(ax2c, [0 log10(clim_val * 2 + 1)]);
        title(ax2c, sprintf('Max accumulator (1-%d)', f), 'Color', 'w');
        xlabel(ax2c, 'Lat [mm]');
        axis(ax2c, 'image');

        sgtitle(fig2, sprintf('SVD cutoff = %d', svdCut), 'Color', 'w', 'FontSize', 16);
        drawnow;
        writeVideo(vw2, getframe(fig2));
        if mod(f, 100) == 0
            fprintf('  %d / %d (%.0f%%)\n', f, nF, f/nF*100);
        end
    end

    close(vw2); close(fig2);
    fprintf('  Saved: %s\n', vidFile2);

    % --- Accumulator snapshot ---
    fig3 = figure('Position', [100 100 1000 800], 'Color', 'k');
    subplot(1,2,1);
    imagesc(xGrid, zGrid, accumMax);
    colormap(gca, 'hot'); colorbar; caxis([0 clim_val * 2]);
    title(sprintf('Max accumulator (linear) - SVD cut=%d, %d frames', svdCut, nF), 'Color', 'w');
    xlabel('Lateral [mm]'); ylabel('Axial [mm]'); axis image;

    subplot(1,2,2);
    imagesc(xGrid, zGrid, log10(accumMax + 1));
    colormap(gca, 'hot'); colorbar;
    title('Max accumulator (log_{10})', 'Color', 'w');
    xlabel('Lateral [mm]'); ylabel('Axial [mm]'); axis image;

    snapFile = fullfile(outputFolder, sprintf('svd_cut%d_accumulator.png', svdCut));
    exportgraphics(fig3, snapFile, 'Resolution', 200, 'BackgroundColor', 'k');
    close(fig3);
    fprintf('\nAccumulator saved: %s\n', snapFile);
end

fprintf('\n=== Done. Outputs in: %s ===\n', outputFolder);
