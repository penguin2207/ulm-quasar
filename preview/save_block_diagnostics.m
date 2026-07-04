function save_block_diagnostics(outDir, blockName, blockIdx, blockResult, ...
    config, xGrid, zGrid, fs_MHz, c, frameRate, pipelineType, varargin)
% SAVE_BLOCK_DIAGNOSTICS  Per-block diagnostic outputs for LAT-ULM/QUASAR.
%
% Saves a multi-panel diagnostic PNG and a summary struct for each block.
% Helps identify noise sources, empty blocks, and per-block quality.
%
% Inputs:
%   outDir       - Output directory (block subfolder created automatically)
%   blockName    - Block filename (used for titles and subfolder)
%   blockIdx     - Block index in the run
%   blockResult  - Struct from process_single_block (LAT-ULM) or per-block
%                  QUASAR data (struct with fields below)
%   config       - Pipeline config struct
%   xGrid, zGrid - Beamforming grid vectors
%   fs_MHz, c    - Sampling freq and speed of sound
%   frameRate    - Compound frame rate [Hz]
%   pipelineType - 'latulm' or 'quasar'
%   varargin     - For QUASAR: {sushiDensity, quasarDensity, velMap, ...
%                                nEnsembles, srX, srZ, svdSpectra}
%
% Outputs: saves to outDir/blocks/block_NNN/
%   block_diagnostic.png  - Multi-panel figure
%   block_summary.mat     - Summary stats

% Create block output directory
shortName = blockName;
if numel(shortName) > 30, shortName = shortName(end-29:end); end
blockDir = fullfile(outDir, 'blocks', sprintf('block_%03d_%s', blockIdx, shortName));
if ~exist(blockDir, 'dir'), mkdir(blockDir); end

switch lower(pipelineType)
    case 'latulm'
        save_latulm_diagnostics(blockDir, blockName, blockIdx, blockResult, ...
            config, xGrid, zGrid, frameRate);
    case 'quasar'
        save_quasar_diagnostics(blockDir, blockName, blockIdx, ...
            config, xGrid, zGrid, frameRate, varargin{:});
end

end

%% ========================================================================
function save_latulm_diagnostics(blockDir, blockName, blockIdx, br, config, xGrid, zGrid, frameRate)

    summary.blockIdx = blockIdx;
    summary.blockName = blockName;
    summary.numFrames = br.numFrames;
    summary.numLocs = br.numLocs;
    summary.numTracks = br.numTracks;
    
    fig = figure('Visible', 'off', 'Position', [0 0 1800 1000]);
    
    % --- Panel 1: Detections per frame over time ---
    subplot(3, 3, 1);
    if ~isempty(br.localizations)
        frames = br.localizations(:, 4);
        [counts, edges] = histcounts(frames, max(1, floor(br.numFrames/50)));
        binCenters = edges(1:end-1) + diff(edges)/2;
        timeMs = binCenters / frameRate * 1000;
        bar(timeMs, counts / diff(edges(1:2)), 'FaceColor', [0.3 0.5 0.8], 'EdgeColor', 'none');
        xlabel('Time [ms]'); ylabel('Detections/frame');
        title(sprintf('Detection Rate\nmean=%.1f/frame', br.numLocs/br.numFrames));
        summary.meanLocsPerFrame = br.numLocs / br.numFrames;
    else
        text(0.5, 0.5, 'No detections', 'HorizontalAlignment', 'center'); axis off;
        summary.meanLocsPerFrame = 0;
    end
    
    % --- Panel 2: Detection spatial distribution (2D histogram) ---
    subplot(3, 3, 2);
    if ~isempty(br.localizations)
        px = 0.050;  % 50 um bins for diagnostics
        xE = min(xGrid):px:max(xGrid)+px;
        zE = min(zGrid):px:max(zGrid)+px;
        dens = histcounts2(br.localizations(:,2), br.localizations(:,1), zE, xE);
        imagesc(xE(1:end-1), zE(1:end-1), log10(dens+1));
        axis image; colormap(gca, hot); colorbar;
        xlabel('Lat [mm]'); ylabel('Ax [mm]');
        title(sprintf('Localization Density\n%d locs', br.numLocs));
    else
        text(0.5, 0.5, 'No localizations', 'HorizontalAlignment', 'center'); axis off;
    end
    
    % --- Panel 3: Detection amplitude histogram ---
    subplot(3, 3, 3);
    if ~isempty(br.localizations)
        amps = br.localizations(:, 3);
        histogram(amps, 50, 'FaceColor', [0.6 0.3 0.3]);
        xlabel('Detection Amplitude'); ylabel('Count');
        title(sprintf('Amplitude Distribution\nmed=%.1f, std=%.1f', median(amps), std(amps)));
        summary.medianAmplitude = median(amps);
        summary.stdAmplitude = std(amps);
    else
        text(0.5, 0.5, 'No data', 'HorizontalAlignment', 'center'); axis off;
        summary.medianAmplitude = 0;
        summary.stdAmplitude = 0;
    end
    
    % --- Panel 4: Tracks ---
    subplot(3, 3, 4);
    if br.numTracks > 0
        hold on;
        cmap = lines(min(br.numTracks, 64));
        for iT = 1:br.numTracks
            t = br.tracks{iT};
            col = cmap(mod(iT-1, size(cmap,1))+1, :);
            plot(t(:,1), t(:,2), '-', 'Color', [col 0.6], 'LineWidth', 0.8);
        end
        set(gca, 'YDir', 'reverse'); axis equal; grid on;
        xlabel('Lat [mm]'); ylabel('Ax [mm]');
        title(sprintf('Tracks (%d)', br.numTracks));
    else
        text(0.5, 0.5, 'No tracks', 'HorizontalAlignment', 'center'); axis off;
    end
    
    % --- Panel 5: Track length histogram ---
    subplot(3, 3, 5);
    if br.numTracks > 0
        lengths = cellfun(@(t) size(t,1), br.tracks);
        histogram(lengths, 20, 'FaceColor', [0.3 0.6 0.3]);
        xlabel('Track Length [frames]'); ylabel('Count');
        title(sprintf('Track Lengths\nmed=%d, max=%d', median(lengths), max(lengths)));
        summary.medianTrackLength = median(lengths);
        summary.maxTrackLength = max(lengths);
    else
        text(0.5, 0.5, 'No tracks', 'HorizontalAlignment', 'center'); axis off;
        summary.medianTrackLength = 0;
        summary.maxTrackLength = 0;
    end
    
    % --- Panel 6: Speed histogram ---
    subplot(3, 3, 6);
    speeds = [];
    if br.numTracks > 0
        for iT = 1:br.numTracks
            t = br.tracks{iT};
            if size(t,1) < 2, continue; end
            dx = diff(t(:,1)); dz = diff(t(:,2)); df = diff(t(:,4));
            dt = df / frameRate;
            sp = sqrt(dx.^2 + dz.^2) ./ (dt + eps);
            speeds = [speeds; sp]; %#ok
        end
    end
    if ~isempty(speeds)
        histogram(speeds, 30, 'FaceColor', [0.8 0.4 0.3]);
        xlim([0 30]); xlabel('Speed [mm/s]'); ylabel('Count');
        title(sprintf('Speed (med=%.1f mm/s)', median(speeds)));
        summary.medianSpeed = median(speeds);
    else
        text(0.5, 0.5, 'No speed data', 'HorizontalAlignment', 'center'); axis off;
        summary.medianSpeed = 0;
    end
    
    % --- Panel 7: Detection threshold sensitivity ---
    subplot(3, 3, 7);
    if isfield(br, 'thresholdCurve') && ~isempty(br.thresholdCurve)
        plot(br.thresholdCurve(:,1), br.thresholdCurve(:,2), 'b.-', 'LineWidth', 1.5);
        xlabel('Threshold'); ylabel('Locs/frame');
        title('Threshold Sensitivity');
        set(gca, 'YScale', 'log'); grid on;
    else
        text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center'); axis off;
    end
    
    % --- Panel 8: Track vs noise classification ---
    subplot(3, 3, 8);
    if ~isempty(br.localizations) && br.numTracks > 0
        % Count localizations that are part of tracks vs orphaned
        trackFrameSets = cellfun(@(t) t(:,4), br.tracks, 'UniformOutput', false);
        trackLocs = 0;
        for iT = 1:br.numTracks
            trackLocs = trackLocs + size(br.tracks{iT}, 1);
        end
        orphanLocs = br.numLocs - trackLocs;
        pie([trackLocs, orphanLocs], {sprintf('Tracked\n%d (%.0f%%)', trackLocs, 100*trackLocs/br.numLocs), ...
            sprintf('Orphaned\n%d (%.0f%%)', orphanLocs, 100*orphanLocs/br.numLocs)});
        title('Localization Fate');
        summary.trackedLocsPercent = 100 * trackLocs / br.numLocs;
    else
        text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center'); axis off;
        summary.trackedLocsPercent = 0;
    end
    
    % --- Panel 9: Per-frame detection count time series ---
    subplot(3, 3, 9);
    if ~isempty(br.localizations)
        allFrames = 1:br.numFrames;
        locsPerFrame = histcounts(br.localizations(:,4), [allFrames, br.numFrames+1]);
        plot(allFrames / frameRate * 1000, locsPerFrame, '.', 'MarkerSize', 2, 'Color', [0.5 0.5 0.5]);
        hold on;
        % Smoothed trend
        if numel(locsPerFrame) > 50
            smoothed = movmean(locsPerFrame, 100);
            plot(allFrames / frameRate * 1000, smoothed, 'r-', 'LineWidth', 1.5);
        end
        xlabel('Time [ms]'); ylabel('Locs/frame');
        title('Detection Time Course');
        ylim([0 max(locsPerFrame)*1.1 + 1]);
        summary.locsStdOverTime = std(locsPerFrame);
    else
        text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center'); axis off;
        summary.locsStdOverTime = 0;
    end
    
    sgtitle(sprintf('Block %d: %s\n%d frames, %d locs, %d tracks', ...
        blockIdx, shortName(blockName), br.numFrames, br.numLocs, br.numTracks), ...
        'FontSize', 11);
    
    saveas(fig, fullfile(blockDir, 'block_diagnostic.png'));
    close(fig);
    
    save(fullfile(blockDir, 'block_summary.mat'), 'summary');
    
    fprintf('      Diagnostics -> %s\n', blockDir);
end

%% ========================================================================
function save_quasar_diagnostics(blockDir, blockName, blockIdx, ...
    config, xGrid, zGrid, frameRate, ...
    blockSushi, blockQuasar, blockVel, nBlockEns, srX, srZ, svdInfo)

    summary.blockIdx = blockIdx;
    summary.blockName = blockName;
    summary.nEnsembles = nBlockEns;
    
    fig = figure('Visible', 'off', 'Position', [0 0 1600 800]);
    
    % --- Panel 1: SUSHI density for this block ---
    subplot(2, 3, 1);
    if any(blockSushi(:) > 0)
        imagesc(srX, srZ, log10(blockSushi / nBlockEns + 1));
        axis image; colormap(gca, hot); colorbar;
        title(sprintf('SUSHI (L1)\n%d ens', nBlockEns));
    else
        text(0.5, 0.5, 'No SUSHI signal', 'HorizontalAlignment', 'center'); axis off;
    end
    xlabel('Lat [mm]'); ylabel('Ax [mm]');
    
    % --- Panel 2: QUASAR density for this block ---
    subplot(2, 3, 2);
    if any(blockQuasar(:) > 0)
        imagesc(srX, srZ, log10(blockQuasar / nBlockEns + 1));
        axis image; colormap(gca, hot); colorbar;
        title('QUASAR (L1+LS)');
    else
        text(0.5, 0.5, 'No QUASAR signal', 'HorizontalAlignment', 'center'); axis off;
    end
    xlabel('Lat [mm]'); ylabel('Ax [mm]');
    
    % --- Panel 3: Velocity for this block ---
    subplot(2, 3, 3);
    if any(blockVel(:) > 0)
        imagesc(srX, srZ, blockVel);
        axis image; colorbar; caxis([0 30]); colormap(gca, jet);
        title('Velocity [mm/s]');
    else
        text(0.5, 0.5, 'No velocity', 'HorizontalAlignment', 'center'); axis off;
    end
    xlabel('Lat [mm]'); ylabel('Ax [mm]');
    
    % --- Panel 4: QUASAR/SUSHI ratio ---
    subplot(2, 3, 4);
    if any(blockSushi(:) > 0) && any(blockQuasar(:) > 0)
        sigMask = blockSushi > 0.01 * max(blockSushi(:));
        rMap = blockQuasar ./ (blockSushi + eps);
        rMap(~sigMask) = NaN;
        imagesc(srX, srZ, rMap); axis image; colorbar; caxis([0 5]);
        colormap(gca, jet);
        validR = rMap(sigMask & ~isnan(rMap));
        title(sprintf('Ratio (mean=%.2f)', mean(validR)));
        summary.meanAmpRatio = mean(validR);
    else
        text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center'); axis off;
        summary.meanAmpRatio = NaN;
    end
    xlabel('Lat [mm]'); ylabel('Ax [mm]');
    
    % --- Panel 5: Lateral profile ---
    subplot(2, 3, 5);
    if any(blockQuasar(:) > 0)
        dP = blockQuasar / nBlockEns;
    elseif any(blockSushi(:) > 0)
        dP = blockSushi / nBlockEns;
    else
        dP = [];
    end
    if ~isempty(dP)
        midZ = round(size(dP,1)/2);
        zBand = max(1,midZ-20):min(size(dP,1),midZ+20);
        lp = mean(dP(zBand,:), 1);
        lp = lp / (max(lp) + eps);
        plot(srX, lp, 'b-', 'LineWidth', 1.5);
        yline(0.5, 'r--', 'FWHM'); grid on;
        xlabel('Lat [mm]'); ylabel('Normalized');
        title('Lateral Profile (mid-depth)');
    else
        text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center'); axis off;
    end
    
    % --- Panel 6: SVD spectrum ---
    subplot(2, 3, 6);
    if ~isempty(svdInfo) && isfield(svdInfo, 'singularValues')
        sv = svdInfo.singularValues;
        semilogy(sv(1:min(50,end)) / sv(1), 'b.-', 'LineWidth', 1.5);
        xlabel('Singular Value Index'); ylabel('Normalized SV');
        title(sprintf('SVD Spectrum\nSV1/SV2=%.1f', sv(1)/sv(2)));
        grid on;
        summary.sv1_sv2_ratio = sv(1) / sv(2);
    else
        text(0.5, 0.5, 'N/A', 'HorizontalAlignment', 'center'); axis off;
        summary.sv1_sv2_ratio = NaN;
    end
    
    % SUSHI support stats
    summary.sushiSupport = sum(blockSushi(:) > 0);
    summary.sushiSupportPct = 100 * summary.sushiSupport / numel(blockSushi);
    summary.sushiMax = max(blockSushi(:));
    summary.quasarMax = max(blockQuasar(:));
    
    sgtitle(sprintf('Block %d: %s\n%d ensembles', blockIdx, shortName(blockName), nBlockEns), ...
        'FontSize', 11);
    
    saveas(fig, fullfile(blockDir, 'block_diagnostic.png'));
    close(fig);
    
    save(fullfile(blockDir, 'block_summary.mat'), 'summary');
    
    fprintf('      Diagnostics -> %s\n', blockDir);
end

%% ========================================================================
function sn = shortName(name)
    if numel(name) > 30, sn = ['...' name(end-27:end)]; else, sn = name; end
end
