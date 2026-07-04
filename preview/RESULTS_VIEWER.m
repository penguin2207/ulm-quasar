%% RESULTS_VIEWER.m
% Comprehensive post-processing viewer for LAT-ULM and QUASAR results.
%
% Features:
%   1. Load any results.mat (auto-detects LAT-ULM vs QUASAR)
%   2. Re-render density/velocity at custom pixel size and zoom region
%   3. Multi-run comparison (concentration series, side-by-side)
%   4. Track statistics (length, speed, direction histograms)
%   5. Channel profile analysis (line cuts, FWHM measurement)
%   6. QUASAR vs SUSHI amplitude comparison
%   7. Export publication-quality figures
%
% Author: Eli Wirth-Apley / Sun Lab, Northeastern University, March 2026

clearvars; close all; clc;

%% ========================================================================
%  CONFIGURATION
%  ========================================================================

% --- Result files to load (add as many as needed) ---
resultFiles = {
    'C:\path\to\ULM3\Primary\LAT_ULM_Results\Batch_5e5\LAT_ULM_results.mat'
    % 'C:\path\to\ULM3\SUSHI-LOW\QUASAR_Results\QUASAR_results.mat'
};

% --- Labels for each file (used in legends and titles) ---
resultLabels = {
    '5x10^5 LAT-ULM'
    % '5x10^5 QUASAR'
};

% --- Display settings ---
display.xLim        = [];       % Lateral zoom [mm]. [] = auto from data
display.zLim        = [];       % Axial zoom [mm]. [] = auto from data
display.srPixel_um  = 50;       % Re-render pixel size [um]
display.velMax      = 30;       % Velocity colorbar max [mm/s]
display.dBrange     = 40;       % Dynamic range for density [dB]
display.exportDPI   = 300;
display.exportDir   = '';       % '' = same as first result file
display.exportFmt   = 'png';

% --- Velocity rendering (track-painted, much cleaner than per-point) ---
display.velSmooth      = true;   % Gaussian smooth velocity map
display.velKernelSize  = 5;      % Kernel size [px]
display.velKernelSigma = 1.0;    % Kernel sigma [px]

% --- Profile analysis ---
profile.enable      = true;
profile.axis        = 'lateral';
profile.position_mm = 17.0;
profile.width_mm    = 0.5;

%% ========================================================================
%  LOAD ALL RESULTS
%  ========================================================================
fprintf('=== RESULTS VIEWER ===\n\n');

nFiles = numel(resultFiles);
data = cell(nFiles, 1);
types = cell(nFiles, 1);

for i = 1:nFiles
    fprintf('[%d/%d] Loading: %s\n', i, nFiles, resultFiles{i});
    if ~exist(resultFiles{i}, 'file')
        fprintf('  WARNING: File not found. Skipping.\n');
        continue;
    end
    
    tmp = load(resultFiles{i});
    
    if isfield(tmp, 'results')
        data{i} = tmp.results;
    elseif isfield(tmp, 'R')
        data{i} = tmp.R;
    else
        fn = fieldnames(tmp);
        data{i} = tmp.(fn{1});
    end
    
    % Auto-detect type
    if isfield(data{i}, 'tracks') || isfield(data{i}, 'localizations')
        types{i} = 'latulm';
    elseif isfield(data{i}, 'sushiDensity') || isfield(data{i}, 'quasarDensity')
        types{i} = 'quasar';
    else
        types{i} = 'unknown';
    end
    fprintf('  Type: %s\n', types{i});
end

% Remove empty entries
valid = ~cellfun(@isempty, data);
data = data(valid); types = types(valid);
resultLabels = resultLabels(valid); resultFiles = resultFiles(valid);
nFiles = numel(data);

if nFiles == 0, error('No valid results files loaded.'); end

% Auto display limits
if isempty(display.xLim) || isempty(display.zLim)
    d = data{1};
    if strcmp(types{1}, 'latulm') && isfield(d, 'tracks') && ~isempty(d.tracks)
        allPts = cell2mat(d.tracks(:));
        if isempty(display.xLim)
            display.xLim = [min(allPts(:,1))-0.5, max(allPts(:,1))+0.5];
        end
        if isempty(display.zLim)
            display.zLim = [min(allPts(:,2))-0.5, max(allPts(:,2))+0.5];
        end
    elseif isfield(d, 'srX')
        if isempty(display.xLim), display.xLim = [d.srX(1), d.srX(end)]; end
        if isempty(display.zLim), display.zLim = [d.srZ(1), d.srZ(end)]; end
    end
end

if isempty(display.exportDir)
    display.exportDir = fileparts(resultFiles{1});
end

latulmIdx = find(strcmp(types, 'latulm'));
quasarIdx = find(strcmp(types, 'quasar'));

%% ========================================================================
%  SECTION 1: DENSITY MAPS
%  ========================================================================
fprintf('\n[1/6] Density maps...\n');

figure('Position', [50 100 400*nFiles 500], 'Color', 'w');

for i = 1:nFiles
    subplot(1, nFiles, i);
    d = data{i};
    
    px = display.srPixel_um / 1000;
    xEdges = display.xLim(1):px:display.xLim(2);
    zEdges = display.zLim(1):px:display.zLim(2);
    xC = xEdges(1:end-1)+px/2;
    zC = zEdges(1:end-1)+px/2;
    
    if strcmp(types{i}, 'latulm')
        density = zeros(numel(zEdges)-1, numel(xEdges)-1);
        tracks = d.tracks;
        for iT = 1:numel(tracks)
            t = tracks{iT};
            for k = 1:size(t,1)
                xi = find(xEdges(1:end-1) <= t(k,1), 1, 'last');
                zi = find(zEdges(1:end-1) <= t(k,2), 1, 'last');
                if ~isempty(xi) && ~isempty(zi) && xi<=size(density,2) && zi<=size(density,1)
                    density(zi,xi) = density(zi,xi) + 1;
                end
            end
        end
        imagesc(xC, zC, log10(density+1));
    elseif strcmp(types{i}, 'quasar') && isfield(d, 'sushiDensity')
        imagesc(d.srX, d.srZ, log10(d.sushiDensity+1));
    end
    
    axis image; colormap(gca, hot); colorbar;
    xlabel('Lateral [mm]'); ylabel('Axial [mm]');
    nLocs=0; nTrk=0; nFr=0;
    if isfield(d,'localizations'), nLocs=size(d.localizations,1); end
    if isfield(d,'tracks'), nTrk=numel(d.tracks); end
    if isfield(d,'totalFrames'), nFr=d.totalFrames; end
    title(sprintf('%s\nDensity: %d locs, %d tracks, %d frames', ...
        resultLabels{i}, nLocs, nTrk, nFr), 'FontSize', 11);
    xlim(display.xLim); ylim(display.zLim);
end

sgtitle('Density Maps');
export_fig_local('density_comparison', display);

%% ========================================================================
%  SECTION 2: VELOCITY MAPS (track-painted, Gaussian-smoothed)
%  ========================================================================
if ~isempty(latulmIdx)
    fprintf('[2/6] Velocity maps (track-painted)...\n');
    
    figure('Position', [50 100 400*numel(latulmIdx) 500], 'Color', 'w');
    
    for ii = 1:numel(latulmIdx)
        i = latulmIdx(ii);
        subplot(1, numel(latulmIdx), ii);
        
        tracks = data{i}.tracks;
        if isempty(tracks), continue; end
        
        frameRate = 1100;
        if isfield(data{i}, 'frameRate'), frameRate = data{i}.frameRate; end
        
        px = display.srPixel_um / 1000;
        xEdges = display.xLim(1):px:display.xLim(2);
        zEdges = display.zLim(1):px:display.zLim(2);
        nXsr = numel(xEdges)-1;
        nZsr = numel(zEdges)-1;
        
        % --- Track-painted velocity ---
        % Paint each track's full trajectory with its median speed.
        % This is much denser and cleaner than per-point speed averaging.
        velAccum = zeros(nZsr, nXsr);
        velWeight = zeros(nZsr, nXsr);
        
        for iT = 1:numel(tracks)
            t = tracks{iT};
            if size(t,1) < 2, continue; end
            
            dx = diff(t(:,1)); dz = diff(t(:,2)); df = diff(t(:,4));
            dt = df / frameRate;
            sp = sqrt(dx.^2 + dz.^2) ./ (dt + eps);
            medSpeed = median(sp);
            
            for k = 1:size(t,1)
                xi = find(xEdges(1:end-1) <= t(k,1), 1, 'last');
                zi = find(zEdges(1:end-1) <= t(k,2), 1, 'last');
                if ~isempty(xi) && ~isempty(zi) && xi<=nXsr && zi<=nZsr
                    velAccum(zi,xi)  = velAccum(zi,xi) + medSpeed;
                    velWeight(zi,xi) = velWeight(zi,xi) + 1;
                end
            end
        end
        
        velMap = velAccum ./ max(velWeight, 1);
        velMap(velWeight == 0) = NaN;
        
        % Optional Gaussian smoothing
        if display.velSmooth
            ks = display.velKernelSize;
            sig = display.velKernelSigma;
            kern = fspecial('gaussian', [ks ks], sig);
            vZ = velMap; vZ(isnan(velMap)) = 0;
            wt = double(~isnan(velMap));
            velMap = conv2(vZ, kern, 'same') ./ max(conv2(wt, kern, 'same'), eps);
            velMap(velWeight == 0) = NaN;
        end
        
        set(gca, 'Color', 'k');  % Black background for empty regions
        hImg = imagesc(xEdges(1:end-1)+px/2, zEdges(1:end-1)+px/2, velMap);
        set(hImg, 'AlphaData', ~isnan(velMap));  % NaN = transparent
        axis image; colorbar; caxis([0 display.velMax]);
        colormap(gca, jet);
        xlabel('Lateral [mm]'); ylabel('Axial [mm]');
        title(sprintf('%s\nVelocity [mm/s]', resultLabels{i}), 'FontSize', 11);
        xlim(display.xLim); ylim(display.zLim);
    end
    sgtitle('Velocity Maps');
    export_fig_local('velocity_comparison', display);
else
    fprintf('[2/6] Velocity maps: no LAT-ULM results, skipping.\n');
end

%% ========================================================================
%  SECTION 3: TRACK STATISTICS
%  ========================================================================
if ~isempty(latulmIdx)
    fprintf('[3/6] Track statistics...\n');
    
    figure('Position', [50 100 1200 400], 'Color', 'w');
    
    for ii = 1:numel(latulmIdx)
        i = latulmIdx(ii);
        tracks = data{i}.tracks;
        if isempty(tracks), continue; end
        
        frameRate = 1100;
        if isfield(data{i}, 'frameRate'), frameRate = data{i}.frameRate; end
        
        lengths = cellfun(@(t) size(t,1), tracks);
        speeds = [];
        for iT = 1:numel(tracks)
            t = tracks{iT};
            if size(t,1) < 2, continue; end
            dx = diff(t(:,1)); dz = diff(t(:,2)); df = diff(t(:,4));
            dt = df / frameRate;
            sp = sqrt(dx.^2 + dz.^2) ./ (dt + eps);
            speeds = [speeds; sp]; %#ok
        end
        
        subplot(1,3,1);
        histogram(lengths, 30, 'FaceColor', [0.4 0.6 0.9]);
        xlabel('Track length [frames]'); ylabel('Count');
        title(sprintf('Track lengths (med=%d)', median(lengths)));
        
        subplot(1,3,2);
        histogram(speeds, 50, 'FaceColor', [0.8 0.4 0.3]);
        xlabel('Speed [mm/s]'); ylabel('Count');
        xline(median(speeds), 'b--', 'LineWidth', 1.5);
        title(sprintf('Speed dist. (med=%.1f mm/s)', median(speeds)));
        
        subplot(1,3,3);
        netDz = cellfun(@(t) t(end,2)-t(1,2), tracks);
        histogram(netDz, 30, 'FaceColor', [0.5 0.8 0.5]);
        xlabel('Net axial displacement [mm]'); ylabel('Count');
        title('Flow direction (>0 = down)');
    end
    
    sgtitle('Track Statistics');
    export_fig_local('track_statistics', display);
else
    fprintf('[3/6] Track statistics: skipping.\n');
end

%% ========================================================================
%  SECTION 4: CHANNEL PROFILE ANALYSIS (FWHM)
%  ========================================================================
if profile.enable
    fprintf('[4/6] Channel profiles...\n');
    
    figure('Position', [100 100 500*nFiles 400], 'Color', 'w');
    
    for i = 1:nFiles
        subplot(1, nFiles, i);
        d = data{i};
        
        px = display.srPixel_um / 1000;
        xEdges = display.xLim(1):px:display.xLim(2);
        zEdges = display.zLim(1):px:display.zLim(2);
        xC = xEdges(1:end-1)+px/2;
        zC = zEdges(1:end-1)+px/2;
        
        if strcmp(types{i}, 'latulm')
            density = zeros(numel(zEdges)-1, numel(xEdges)-1);
            tracks = d.tracks;
            for iT = 1:numel(tracks)
                t = tracks{iT};
                for k = 1:size(t,1)
                    xi = find(xEdges(1:end-1) <= t(k,1), 1, 'last');
                    zi = find(zEdges(1:end-1) <= t(k,2), 1, 'last');
                    if ~isempty(xi) && ~isempty(zi) && xi<=size(density,2) && zi<=size(density,1)
                        density(zi,xi) = density(zi,xi) + 1;
                    end
                end
            end
            posAxis = xC;
            if strcmp(profile.axis, 'lateral')
                zMask = abs(zC - profile.position_mm) <= profile.width_mm/2;
                profileLine = mean(density(zMask, :), 1);
            else
                xMask = abs(xC - profile.position_mm) <= profile.width_mm/2;
                profileLine = mean(density(:, xMask), 2);
                posAxis = zC;
            end
        elseif strcmp(types{i}, 'quasar') && isfield(d, 'sushiDensity')
            densityMap = d.sushiDensity;
            if strcmp(profile.axis, 'lateral')
                zMask = abs(d.srZ - profile.position_mm) <= profile.width_mm/2;
                profileLine = mean(densityMap(zMask, :), 1);
                posAxis = d.srX;
            else
                xMask = abs(d.srX - profile.position_mm) <= profile.width_mm/2;
                profileLine = mean(densityMap(:, xMask), 2);
                posAxis = d.srZ;
            end
        else
            continue;
        end
        
        profileLine = profileLine / (max(profileLine) + eps);
        plot(posAxis, profileLine, 'b-', 'LineWidth', 1.5);
        hold on;
        yline(0.5, 'r--', 'FWHM', 'LineWidth', 1);
        xlabel(sprintf('%s [mm]', profile.axis)); ylabel('Normalized');
        title(sprintf('%s\n%s @ %.1f mm', resultLabels{i}, ...
            profile.axis, profile.position_mm), 'FontSize', 10);
        grid on;
        
        % Measure FWHM for each peak
        aboveHalf = profileLine >= 0.5;
        regions = diff([0, aboveHalf, 0]);
        starts = find(regions == 1);
        ends = find(regions == -1) - 1;
        
        for r = 1:numel(starts)
            peakRegion = posAxis(starts(r):ends(r));
            fwhmVal = peakRegion(end) - peakRegion(1);
            peakCenter = mean(peakRegion);
            text(peakCenter, 0.55, sprintf('%.0f um', fwhmVal*1000), ...
                'HorizontalAlignment', 'center', 'FontSize', 9, ...
                'Color', 'r', 'FontWeight', 'bold');
        end
    end
    
    sgtitle(sprintf('Channel Profiles (%s at %.1f mm)', profile.axis, profile.position_mm));
    export_fig_local('channel_profiles', display);
else
    fprintf('[4/6] Channel profiles: disabled.\n');
end

%% ========================================================================
%  SECTION 5: QUASAR-SPECIFIC VIEWS
%  ========================================================================
if ~isempty(quasarIdx)
    fprintf('[5/6] QUASAR views...\n');
    for ii = 1:numel(quasarIdx)
        i = quasarIdx(ii);
        d = data{i};
        if isfield(d, 'sushiDensity') && isfield(d, 'quasarDensity')
            figure('Position', [100 100 1200 400], 'Color', 'w');
            
            subplot(1,3,1);
            imagesc(d.srX, d.srZ, log10(d.sushiDensity+1));
            axis image; colormap(gca, hot); colorbar;
            title('SUSHI (L1 only)'); xlabel('Lateral [mm]'); ylabel('Axial [mm]');
            
            subplot(1,3,2);
            imagesc(d.srX, d.srZ, log10(d.quasarDensity+1));
            axis image; colormap(gca, hot); colorbar;
            title('QUASAR (L1 + LS)'); xlabel('Lateral [mm]'); ylabel('Axial [mm]');
            
            subplot(1,3,3);
            ratio = d.quasarDensity ./ max(d.sushiDensity, eps);
            imagesc(d.srX, d.srZ, ratio);
            axis image; colormap(gca, jet); colorbar; caxis([0 5]);
            title('QUASAR / SUSHI ratio'); xlabel('Lateral [mm]'); ylabel('Axial [mm]');
            
            sgtitle(resultLabels{i});
            export_fig_local(sprintf('quasar_comparison_%d', ii), display);
        end
    end
else
    fprintf('[5/6] QUASAR views: no QUASAR results.\n');
end

%% ========================================================================
%  SECTION 6: MULTI-RUN OVERLAY
%  ========================================================================
if nFiles > 1
    fprintf('[6/6] Multi-run overlay...\n');
    fprintf('  (See density_comparison figure)\n');
else
    fprintf('[6/6] Multi-run overlay: only 1 file, skipping.\n');
end

fprintf('\n=== RESULTS VIEWER COMPLETE ===\n');
fprintf('Figures saved to: %s\n', display.exportDir);

%% ========================================================================
%  HELPER
%  ========================================================================
function export_fig_local(name, display)
    fname = fullfile(display.exportDir, [name '.' display.exportFmt]);
    try
        exportgraphics(gcf, fname, 'Resolution', display.exportDPI);
    catch
        saveas(gcf, fname);
    end
    if strcmp(display.exportFmt, 'png')
        try savefig(gcf, fullfile(display.exportDir, [name '.fig'])); catch, end
    end
    fprintf('  Saved: %s\n', fname);
end
