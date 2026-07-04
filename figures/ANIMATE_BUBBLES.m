%% ANIMATE_BUBBLES.m
% Animate microbubble localizations and tracks from LAT-ULM results.
%
% Generates a video showing bubble positions evolving over time, with
% trailing tracks and optional velocity coloring.
%
% Features:
%   - Track smoothing (movmean) to reduce visual jitter
%   - Speed, direction, or per-track coloring
%   - Configurable clip duration (avoids overly long renders)
%   - Static three-panel summary figure
%
% Author: Eli Wirth-Apley / Sun Lab, Northeastern University, March 2026

clearvars; close all; clc;

%% CONFIGURATION
% Apr 17 batch format (default): LAT_ULM/<C>_svd<N>_latulm.mat + batch_config.mat
% Legacy format supported too: LAT_ULM_results.mat with R.localizations/.tracks/.frameRate.
resultsFile = 'E:\4-17-Final\Results\LAT_ULM\C8_svd5_latulm.mat';
branch      = 'combined';       % 'single' (b2 only) or 'combined' (b1+b2+b3). Apr 17 only.
outputVideo = 'E:\4-17-Final\Results\Figures\C8_combined_15sec.mp4';

% Display settings
tailLength    = 60;     % Show this many trailing frames of each track
frameStep     = 30;     % Advance this many frames per animation frame. For 15s output at 30 fps on ~13851 combined frames, use ~30.
markerSize    = 8;      % Bubble marker size
trailAlpha    = 0.4;    % Trail transparency
fps           = 30;     % Output video frame rate
dpi           = 150;    % Output resolution

% Jitter reduction: moving average on track x/z positions
smoothWindow  = 5;      % Odd number. 1 = no smoothing. 5-7 works well.

% Minimum speed filter: remove stationary/wall-stuck tracks
minSpeed_mm_s = 2.0;    % Reject tracks with median speed below this [mm/s]. 0 = no filter.

% Clip duration: set maxDuration_s to limit video length ([] = full dataset)
maxDuration_s = 15;     % Output video max length [s]. [] = no limit.

% Focus region (set to [] for auto from results)
xLim          = [-9, -4];    % Crop tightly to the C8 tube region
zLim          = [14, 22];

% Color mode: 'track' (color per track), 'speed', or 'direction'
colorMode     = 'speed';    % Poiseuille profile shows up as a color spread through the tube

%% Load results
% Detect schema: Apr 17 batch saves `latulmResult.{single,combined}` in
% per-(C, SVD) files plus a sibling `batch_config.mat` for frameRate.
% Legacy format had top-level `R.localizations/.tracks/.frameRate`.
fprintf('Loading results...\n');
tmp = load(resultsFile);

if isfield(tmp, 'latulmResult')
    % Apr 17 schema
    LR = tmp.latulmResult;
    if ~isfield(LR, branch)
        error('Apr 17 result has no branch ''%s''. Use ''single'' or ''combined''.', branch);
    end
    br = LR.(branch);
    R.tracks = br.tracks;
    R.localizations = br.locs;       % [x_mm, z_mm, amplitude, globalFrameIdx]

    % Pull frameRate from the sibling batch_config.mat (Results dir, one up
    % from LAT_ULM/).
    [latulmDir, ~, ~] = fileparts(resultsFile);
    [resultsDir, ~, ~] = fileparts(latulmDir);
    cfgFile = fullfile(resultsDir, 'batch_config.mat');
    if ~exist(cfgFile, 'file')
        error('Apr 17 schema requires batch_config.mat at %s', cfgFile);
    end
    BC = load(cfgFile, 'meta');
    R.frameRate = BC.meta.frameRate_Hz;

    % totalFrames: max frame index seen in any track (b1+b2+b3 concatenation
    % means combined goes up to ~3 * single)
    R.totalFrames = 0;
    for ii = 1:numel(R.tracks)
        if ~isempty(R.tracks{ii})
            R.totalFrames = max(R.totalFrames, max(R.tracks{ii}(:,4)));
        end
    end

    fprintf('  Apr 17 schema: %s/%s, conc=%s (%.1e), SVD=%d\n', ...
        branch, LR.concLabel, LR.concLabel, LR.concValue, LR.svd);

elseif isfield(tmp, 'results')
    R = tmp.results;
elseif isfield(tmp, 'R')
    R = tmp.R;
else
    fn = fieldnames(tmp);
    R = tmp.(fn{1});
end
clear tmp;

locs = R.localizations;      % [x_mm, z_mm, amplitude, globalFrameIdx]
tracks = R.tracks;           % Cell array, each [Nx4]: [x, z, amp, frame]
frameRate = R.frameRate;
totalFrames = 0;
if isfield(R, 'totalFrames'), totalFrames = R.totalFrames; end

fprintf('  %d localizations, %d tracks, %d frames at %.0f Hz\n', ...
    size(locs,1), numel(tracks), totalFrames, frameRate);

if isempty(tracks)
    error('No tracks in results. Nothing to animate.');
end

nTracks = numel(tracks);

%% Filter slow/stationary tracks
if minSpeed_mm_s > 0
    keepMask = true(nTracks, 1);
    for i = 1:nTracks
        t = tracks{i};
        if size(t,1) > 1
            dx = diff(t(:,1)); dz = diff(t(:,2)); df = diff(t(:,4));
            dt = df / frameRate;
            sp = sqrt(dx.^2 + dz.^2) ./ (dt + eps);
            if median(sp) < minSpeed_mm_s
                keepMask(i) = false;
            end
        else
            keepMask(i) = false;
        end
    end
    nRemoved = sum(~keepMask);
    tracks = tracks(keepMask);
    nTracks = numel(tracks);
    fprintf('  Filtered %d slow tracks (< %.1f mm/s), %d remaining\n', ...
        nRemoved, minSpeed_mm_s, nTracks);
end

%% Smooth tracks to reduce jitter
if smoothWindow > 1
    fprintf('Smoothing tracks (window=%d)...\n', smoothWindow);
    for i = 1:nTracks
        t = tracks{i};
        if size(t,1) > smoothWindow
            tracks{i}(:,1) = movmean(t(:,1), smoothWindow);
            tracks{i}(:,2) = movmean(t(:,2), smoothWindow);
            % Keep amplitude and frame index untouched
        end
    end
end

%% Build frame index
fprintf('Building frame index...\n');
trackFrameRanges = zeros(nTracks, 2);
for i = 1:nTracks
    trackFrameRanges(i, :) = [min(tracks{i}(:,4)), max(tracks{i}(:,4))];
end

% Auto-determine display limits from track data
allTrackPts = cell2mat(tracks(:));
if isempty(xLim)
    xLim = [min(allTrackPts(:,1))-0.5, max(allTrackPts(:,1))+0.5];
end
if isempty(zLim)
    zLim = [min(allTrackPts(:,2))-0.5, max(allTrackPts(:,2))+0.5];
end

% Frame range with actual track data
frameMin = min(trackFrameRanges(:,1));
frameMax = max(trackFrameRanges(:,2));

% Clip to max duration
if ~isempty(maxDuration_s)
    maxDataFrames = maxDuration_s * fps * frameStep;
    if (frameMax - frameMin) > maxDataFrames
        frameMax = frameMin + maxDataFrames;
        fprintf('  Clipping to %.0f s of data\n', maxDuration_s);
    end
end

animFrames = frameMin:frameStep:frameMax;
fprintf('  Animating frames %d to %d (step=%d, %d animation frames, %.1f s)\n', ...
    frameMin, frameMax, frameStep, numel(animFrames), numel(animFrames)/fps);

%% Precompute track colors and speeds
% Per-track colors
cmap = lines(min(nTracks, 256));
trackColors = cmap(mod((1:nTracks)-1, size(cmap,1))+1, :);

% Per-track speeds and directions
trackMedSpeeds = zeros(nTracks, 1);
trackNetDir    = zeros(nTracks, 1);
trackSpeeds    = cell(nTracks, 1);

for i = 1:nTracks
    t = tracks{i};
    if size(t,1) > 1
        dx = diff(t(:,1)); dz = diff(t(:,2));
        df = diff(t(:,4));
        dt = df / frameRate;
        sp = sqrt(dx.^2 + dz.^2) ./ (dt + eps);
        trackSpeeds{i}    = [sp(1); sp];
        trackMedSpeeds(i) = median(sp);
        trackNetDir(i)    = sign(t(end,2) - t(1,2));
    else
        trackSpeeds{i} = 0;
    end
end

speedMax = prctile(trackMedSpeeds(trackMedSpeeds > 0), 95);
if speedMax == 0, speedMax = 20; end
fprintf('  Speed range: 0 to %.1f mm/s (95th percentile)\n', speedMax);

%% Create animation
fprintf('Generating animation...\n');

fig = figure('Position', [100 100 800 900], 'Color', 'k', 'MenuBar', 'none');
ax = axes('Parent', fig, 'Color', 'k', 'XColor', 'w', 'YColor', 'w');
hold(ax, 'on');
% Lock FOV: explicit data aspect + manual limit modes prevent the per-frame
% jitter that `axis equal` produces (auto-rebalancing on new plot data).
set(ax, 'YDir', 'reverse', ...
        'DataAspectRatio', [1 1 1], ...
        'PlotBoxAspectRatioMode', 'manual', ...
        'XLim', xLim, 'YLim', zLim, ...
        'XLimMode', 'manual', 'YLimMode', 'manual');
xlabel(ax, 'Lateral [mm]', 'Color', 'w');
ylabel(ax, 'Axial [mm]', 'Color', 'w');
set(ax, 'FontSize', 12);

% Video writer
[vidPath, vidName, vidExt] = fileparts(outputVideo);
if isempty(vidPath), vidPath = pwd; end
vidFile = fullfile(vidPath, [vidName vidExt]);
v = VideoWriter(vidFile, 'MPEG-4');
v.FrameRate = fps;
v.Quality = 95;
open(v);

wb = waitbar(0, 'Rendering animation...', 'Name', 'Bubble Animation');

for iAnim = 1:numel(animFrames)
    currentFrame = animFrames(iAnim);
    
    if mod(iAnim, 50) == 0 && ishandle(wb)
        waitbar(iAnim/numel(animFrames), wb, ...
            sprintf('Frame %d/%d (data frame %d)', iAnim, numel(animFrames), currentFrame));
    end
    
    cla(ax);
    % Re-lock FOV after cla in case any prior plot triggered rescaling
    set(ax, 'XLim', xLim, 'YLim', zLim, ...
            'XLimMode', 'manual', 'YLimMode', 'manual');

    windowStart = currentFrame - tailLength;
    windowEnd = currentFrame;
    
    activeMask = (trackFrameRanges(:,1) <= windowEnd) & (trackFrameRanges(:,2) >= windowStart);
    activeIdx = find(activeMask);
    
    for ii = 1:numel(activeIdx)
        ti = activeIdx(ii);
        t = tracks{ti};
        
        inWindow = (t(:,4) >= windowStart) & (t(:,4) <= windowEnd);
        tw = t(inWindow, :);
        if isempty(tw), continue; end
        
        % Determine color
        switch colorMode
            case 'speed'
                frac = min(trackMedSpeeds(ti) / speedMax, 1);
                col = [frac, 0.1, 1-frac];  % blue -> red
                
            case 'direction'
                if trackNetDir(ti) > 0
                    col = [1 0.3 0.3];    % Red = downward
                elseif trackNetDir(ti) < 0
                    col = [0.3 0.3 1];    % Blue = upward
                else
                    col = [0.5 0.5 0.5];
                end
                
            otherwise  % 'track'
                col = trackColors(ti, :);
        end
        
        % Draw trail
        if size(tw, 1) > 1
            plot(ax, tw(:,1), tw(:,2), '-', 'Color', [col, trailAlpha], 'LineWidth', 1.5);
        end
        
        % Draw current position
        [~, latestIdx] = max(tw(:,4));
        plot(ax, tw(latestIdx,1), tw(latestIdx,2), '.', ...
            'Color', col, 'MarkerSize', markerSize);
    end
    
    % Time annotation
    timeS = (currentFrame - frameMin) / frameRate;
    title(ax, sprintf('Frame %d | t = %.2f s | %d active tracks', ...
        currentFrame, timeS, numel(activeIdx)), 'Color', 'w', 'FontSize', 11);
    
    drawnow limitrate;
    frame = getframe(fig);
    writeVideo(v, frame);
end

close(wb);
close(v);
close(fig);

fprintf('\nAnimation saved to: %s\n', vidFile);
fprintf('  Duration: %.1f s at %d fps\n', numel(animFrames)/fps, fps);
fprintf('  Data time covered: %.1f s\n', (frameMax-frameMin)/frameRate);

%% Static summary figure (three-panel)
figure('Position', [100 100 1200 500], 'Color', 'w');

% Panel 1: All tracks colored by track ID
subplot(1,3,1);
hold on;
for i = 1:nTracks
    t = tracks{i};
    col = trackColors(min(i, size(trackColors,1)), :);
    plot(t(:,1), t(:,2), '-', 'Color', [col 0.6], 'LineWidth', 0.8);
end
set(gca, 'YDir', 'reverse'); axis equal;
xlim(xLim); ylim(zLim);
xlabel('Lateral [mm]'); ylabel('Axial [mm]');
title(sprintf('All tracks (%d)', nTracks));

% Panel 2: Tracks colored by mean speed
subplot(1,3,2);
hold on;
speedLim = prctile(trackMedSpeeds(trackMedSpeeds > 0), 95);
if speedLim == 0, speedLim = 1; end

for i = 1:nTracks
    t = tracks{i};
    frac = min(trackMedSpeeds(i) / speedLim, 1);
    col = [frac, 0.1, 1-frac];
    plot(t(:,1), t(:,2), '-', 'Color', [col 0.6], 'LineWidth', 0.8);
end
set(gca, 'YDir', 'reverse'); axis equal;
xlim(xLim); ylim(zLim);
xlabel('Lateral [mm]'); ylabel('Axial [mm]');
title(sprintf('Tracks by speed (0-%.0f mm/s)', speedLim));
colormap(gca, [linspace(0,1,64)', 0.1*ones(64,1), linspace(1,0,64)']);
cb = colorbar; caxis([0 speedLim]);
ylabel(cb, 'Speed [mm/s]');

% Panel 3: Tracks colored by direction
subplot(1,3,3);
hold on;
for i = 1:nTracks
    t = tracks{i};
    if trackNetDir(i) > 0
        col = [1 0.3 0.3];
    elseif trackNetDir(i) < 0
        col = [0.3 0.3 1];
    else
        col = [0.5 0.5 0.5];
    end
    plot(t(:,1), t(:,2), '-', 'Color', [col 0.6], 'LineWidth', 0.8);
end
set(gca, 'YDir', 'reverse'); axis equal;
xlim(xLim); ylim(zLim);
xlabel('Lateral [mm]'); ylabel('Axial [mm]');
title('Flow direction (red=down, blue=up)');

sgtitle(sprintf('LAT-ULM Track Analysis: %d tracks, %.0f Hz', nTracks, frameRate));
saveas(gcf, fullfile(fileparts(vidFile), 'track_analysis.png'));

fprintf('\nStatic figure saved to: track_analysis.png\n');
