function blockResult = process_single_block(dataFolder, baseFilename, modeName, ...
    config, anglePairs, elemPos_mm, delayTables, xGrid, zGrid, nX, nZ, fs_MHz, c, depthOffset_mm, bgMeanIQ)
% PROCESS_SINGLE_BLOCK  Run full LAT-ULM processing on one VADA block.
%
% Handles: load → beamform (GPU) → PI compound → SVD filter → detect →
%          localize → track. Returns all results for accumulation.
%
% Input:
%   dataFolder, baseFilename, modeName - VADA file identifiers
%   config       - Full pipeline config struct
%   anglePairs   - Angle/polarity mapping from metadata
%   elemPos_mm   - Element positions [mm]
%   delayTables  - Precomputed GPU delay tables (cell array, one per angle)
%   xGrid, zGrid - Beamforming grid axes [mm]
%   nX, nZ       - Grid dimensions
%   fs_MHz, c, depthOffset_mm - Acoustic parameters
%
% Output:
%   blockResult  - Struct with fields:
%     .localizations  - [N x 4]: [x_mm, z_mm, amplitude, frameIdx]
%     .tracks         - Cell array of tracks (each [M x 4])
%     .timestamps     - [nFrames x 1] timestamps [ms]
%     .numFrames      - Number of compound frames processed
%     .numLocs        - Total localizations
%     .numTracks      - Total tracks
%     .processingTime - Elapsed time [s]
%     .filename       - Source filename

tBlock = tic;

% Initialize output
blockResult.filename  = baseFilename;
blockResult.localizations = [];
blockResult.tracks    = {};
blockResult.timestamps = [];
blockResult.numFrames = 0;
blockResult.numLocs   = 0;
blockResult.numTracks = 0;

% --- Determine block size ---
try
    [VadaTest, ~, ~, BlockConfig] = VsiVadaDataRead(dataFolder, baseFilename, ...
        1:config.eventsPerFrame, modeName);
    clear VadaTest;
    numTotalEvents = numel(BlockConfig.PulseSequences(1).Events);
catch ME
    fprintf('    ERROR reading %s: %s\n', baseFilename, ME.message);
    blockResult.processingTime = toc(tBlock);
    return;
end

numCompoundFrames = floor(numTotalEvents / config.eventsPerFrame);
numChunks = ceil(numCompoundFrames / config.chunkSize);

fprintf('    %s: %d events -> %d frames, %d chunks\n', ...
    baseFilename, numTotalEvents, numCompoundFrames, numChunks);

% --- Chunk-wise processing ---
% Accumulate localizations in a cell array (O(1) append), then concat
% once at the end. Avoids the O(n^2) growth of `[A; newRows]` in a loop.
locsPerChunk = cell(numChunks, 1);
allTimestamps = zeros(numCompoundFrames, 1);
tsFillIdx = 0;
globalFrameIdx = 0;

% Phase timers (accumulated across chunks)
t_load = 0; t_bf = 0; t_mc = 0; t_bg = 0; t_svd = 0; t_det = 0;
t_detectFn = 0; t_centroid = 0;  % sub-timers inside the det phase

for iChunk = 1:numChunks
    frameStart = (iChunk-1)*config.chunkSize + 1;
    frameEnd   = min(iChunk*config.chunkSize, numCompoundFrames);
    nFrames    = frameEnd - frameStart + 1;
    eventStart = (frameStart-1)*config.eventsPerFrame + 1;
    eventEnd   = frameEnd * config.eventsPerFrame;

    % Load
    tLoad = tic;
    try
        [VadaChunk,~,~,~] = VsiVadaDataRead(dataFolder, baseFilename, ...
            eventStart:eventEnd, modeName);
    catch ME
        fprintf('      Chunk %d/%d SKIPPED: %s\n', iChunk, numChunks, ME.message);
        continue;
    end
    t_load = t_load + toc(tLoad);
    
    % Beamform + PI + Compound (or zero-only)
    %
    % When useCUDA is on, we BATCH all frames of a chunk into a single
    % beamformer call (one per angle). This amortizes the per-call
    % overhead (cuFFT plan creation + cudaMalloc/Free) across hundreds
    % of frames instead of paying it every frame. Huge speedup for
    % multi-angle compound mode.
    %
    % When useCUDA is off, fall back to per-frame calls because the
    % MATLAB gpuArray beamformer only accepts 2D input.
    tBf = tic;
    timestamps = zeros(nFrames, 1);

    useBatch = isfield(config, 'useCUDA') && config.useCUDA;
    if useBatch
        bf_fn = @beamform_cuda;
    else
        bf_fn = @beamform_planewave_gpu;
    end

    % Collect timestamps up front
    for iFrame = 1:nFrames
        baseEvt = (iFrame-1)*config.eventsPerFrame;
        timestamps(iFrame) = VadaChunk(baseEvt+1).Timestamp;
    end

    if useBatch
        % --- Batched path: one MEX call per angle for the whole chunk ---
        nSamplesRf = size(VadaChunk(1).Data, 1);
        nRxRf      = size(VadaChunk(1).Data, 2);

        if config.zeroOnly
            % Collect all frames for the zero-degree event
            rfBatch = zeros(nSamplesRf, nRxRf, nFrames, 'single');
            for iFrame = 1:nFrames
                baseEvt = (iFrame-1)*config.eventsPerFrame;
                evIdx = baseEvt + config.zeroEventIdx;
                rfBatch(:,:,iFrame) = single(VadaChunk(evIdx).Data);
            end
            rxPos_mm = elemPos_mm(anglePairs(config.zeroAngleIdx).rxElements);
            IQ_compound = bf_fn(rfBatch, rxPos_mm, ...
                anglePairs(config.zeroAngleIdx).angle, ...
                xGrid, zGrid, fs_MHz, c, depthOffset_mm, delayTables{1});
            IQ_compound = single(IQ_compound);
            clear rfBatch;
        else
            % Multi-angle: batch per angle, accumulate complex images
            IQ_compound = complex(zeros(nZ, nX, nFrames, 'single'));
            rfBatch = zeros(nSamplesRf, nRxRf, nFrames, 'single');
            for a = 1:config.numAngles
                for iFrame = 1:nFrames
                    baseEvt = (iFrame-1)*config.eventsPerFrame;
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
                        rfBatch(:,:,iFrame) = rfPos + rfNeg;  % PI summation
                    else
                        rfBatch(:,:,iFrame) = rfPos;
                    end
                end

                rxPos_mm = elemPos_mm(anglePairs(a).rxElements);
                bfBatch = bf_fn(rfBatch, rxPos_mm, anglePairs(a).angle, ...
                    xGrid, zGrid, fs_MHz, c, depthOffset_mm, delayTables{a});
                IQ_compound = IQ_compound + single(bfBatch);
            end
            clear rfBatch bfBatch;
        end
    else
        % --- Per-frame path (MATLAB gpuArray beamformer, 2D only) ---
        IQ_compound = zeros(nZ, nX, nFrames, 'single');
        for iFrame = 1:nFrames
            baseEvt = (iFrame-1)*config.eventsPerFrame;

            if config.zeroOnly
                evIdx = baseEvt + config.zeroEventIdx;
                rfData = single(VadaChunk(evIdx).Data);
                rxPos_mm = elemPos_mm(anglePairs(config.zeroAngleIdx).rxElements);
                bfImg = bf_fn(rfData, rxPos_mm, ...
                    anglePairs(config.zeroAngleIdx).angle, ...
                    xGrid, zGrid, fs_MHz, c, depthOffset_mm, delayTables{1});
                IQ_compound(:,:,iFrame) = single(bfImg);
            else
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
                    rxPos_mm = elemPos_mm(anglePairs(a).rxElements);
                    bfImg = bf_fn(rfData, rxPos_mm, anglePairs(a).angle, ...
                        xGrid, zGrid, fs_MHz, c, depthOffset_mm, delayTables{a});
                    compImg = compImg + single(bfImg);
                end
                IQ_compound(:,:,iFrame) = compImg;
            end
        end
    end
    clear VadaChunk;
    t_bf = t_bf + toc(tBf);

    % Motion correction (if enabled) — must precede SVD for clean clutter separation
    if isfield(config, 'motionCorrection') && config.motionCorrection.enable
        tMc = tic;
        mcParams.method    = config.motionCorrection.method;
        mcParams.refType   = config.motionCorrection.refType;
        mcParams.refWindow = config.motionCorrection.refWindow;
        mcParams.maxShift  = config.motionCorrection.maxShift;

        [mcShifts, mcDiag] = estimate_tissue_motion(IQ_compound, mcParams);
        IQ_compound = apply_motion_correction(IQ_compound, mcShifts, config.useGPU);
        t_mc = t_mc + toc(tMc);

        if iChunk == 1
            fprintf('    Motion correction: max=%.2f px, mean=%.2f px\n', ...
                mcDiag.maxDisp_px, mcDiag.meanDisp_px);
        end
    end

    % Background subtraction (if provided)
    if ~isempty(bgMeanIQ)
        tBg = tic;
        for iFrame = 1:nFrames
            IQ_compound(:,:,iFrame) = IQ_compound(:,:,iFrame) - bgMeanIQ;
        end
        t_bg = t_bg + toc(tBg);
    end

    % SVD clutter filter (randomized SVD if enabled and cutoffHigh not set)
    tSvd = tic;
    if isfield(config, 'useRSVD') && config.useRSVD && ...
       (isempty(config.svd.cutoffHigh) || config.svd.cutoffHigh == 0)
        IQ_filtered = svd_clutter_filter_rsvd(IQ_compound, ...
            config.svd.cutoffLow, config.svd.cutoffHigh, config.useGPU);
    else
        IQ_filtered = svd_clutter_filter_gpu(IQ_compound, ...
            config.svd.cutoffLow, config.svd.cutoffHigh, config.useGPU);
    end
    clear IQ_compound;
    t_svd = t_svd + toc(tSvd);

    % Detection + Localization (vectorized)
    %
    % Per-frame: call detect_microbubbles, then batch ALL detections'
    % centroids into one vectorized computation:
    %   1. Zero-pad the frame so every ROI is full-size (no per-ROI clamping).
    %   2. Build a [roiSize x roiSize x nDets] index array via broadcasting.
    %   3. Extract all ROIs in one indexing op.
    %   4. Flatten to [roiPx x nDets] and compute all centroids as a single
    %      matrix-vector product.
    %
    % chunkLocs is pre-allocated (avoids the O(n^2) dynamic growth that
    % dominated when a chunk produced thousands of localizations).
    tDet = tic;
    halfROI = floor(config.det.roiSize_px / 2);
    roiSize = 2*halfROI + 1;
    roiPx   = roiSize * roiSize;

    % Pre-compute centroid weight vectors (once per chunk)
    [colGrid, rowGrid] = meshgrid(1:roiSize, 1:roiSize);
    rowVec = single(rowGrid(:));   % [roiPx x 1]
    colVec = single(colGrid(:));

    % Pre-allocate chunkLocs with a generous initial capacity
    chunkLocsCap = max(nFrames * 500, 10000);
    chunkLocs    = zeros(chunkLocsCap, 4);
    nChunkLocs   = 0;

    for iFrame = 1:nFrames
        tFn = tic;
        frame = abs(IQ_filtered(:,:,iFrame));
        dets  = detect_microbubbles(frame, config.det, config.bf.dx, config.bf.dz);
        t_detectFn = t_detectFn + toc(tFn);
        if isempty(dets), continue; end
        nDets = size(dets, 1);

        tCen = tic;
        % Zero-pad so every ROI is a full (roiSize x roiSize) block.
        % This lets us extract all ROIs via a single indexing op instead
        % of a per-detection clamp-and-slice loop.
        frame_padded = padarray(single(frame), [halfROI halfROI], 0, 'both');
        nZp = size(frame_padded, 1);

        % Detection centers in padded coordinates
        det_r_pad = int32(dets(:,1)) + int32(halfROI);
        det_c_pad = int32(dets(:,2)) + int32(halfROI);

        % Build [roiSize x roiSize x nDets] linear-index array via implicit
        % expansion (broadcasting):
        %   r_idx(i, 1, k) = det_r_pad(k) + offsets(i)
        %   c_idx(1, j, k) = det_c_pad(k) + offsets(j)
        %   lin   (i, j, k) = r_idx(i,1,k) + (c_idx(1,j,k) - 1) * nZp
        offsets_col = int32(-halfROI:halfROI).';              % [roiSize x 1]
        offsets_row = reshape(offsets_col, 1, roiSize);       % [1 x roiSize]
        r_idx = reshape(det_r_pad, 1, 1, nDets) + offsets_col;      % [roiSize x 1 x nDets]
        c_idx = reshape(det_c_pad, 1, 1, nDets) + offsets_row;      % [1 x roiSize x nDets]

        lin_idx = r_idx + (c_idx - 1) * int32(nZp);           % [roiSize x roiSize x nDets]
        rois    = frame_padded(lin_idx);                       % [roiSize x roiSize x nDets]

        % Flatten to [roiPx x nDets] for matrix-level centroid math
        rois_flat = reshape(rois, roiPx, nDets);

        % Background subtract per ROI (matches intensity_weighted_centroid)
        roi_mins = min(rois_flat, [], 1);                     % [1 x nDets]
        rois_bg  = max(rois_flat - roi_mins, 0);              % [roiPx x nDets]

        % Weighted centroid: two matrix-vector products
        tots  = sum(rois_bg, 1);                              % [1 x nDets]
        sum_r = rowVec' * rois_bg;                            % [1 x nDets]
        sum_c = colVec' * rois_bg;

        % Safe division (fallback to ROI center when total = 0)
        safe_tots = max(tots, eps('single'));
        subR = sum_r ./ safe_tots;
        subC = sum_c ./ safe_tots;
        zero_mask = (tots == 0);
        center_val = single((roiSize + 1) / 2);
        subR(zero_mask) = center_val;
        subC(zero_mask) = center_val;

        % Max intensity from original (non-bg-subtracted) ROI
        maxVals = max(rois_flat, [], 1);                      % [1 x nDets]

        % Convert to mm positions. Use detection center + offset from ROI
        % center; this matches the original formula zGrid(r1)+(sR-1)*dz
        % for interior detections (r1 = det_r - halfROI).
        x_center = reshape(xGrid(dets(:,2)), nDets, 1);
        z_center = reshape(zGrid(dets(:,1)), nDets, 1);
        subC_col = double(subC(:));
        subR_col = double(subR(:));
        x_pos = x_center + (subC_col - halfROI - 1) * config.bf.dx;
        z_pos = z_center + (subR_col - halfROI - 1) * config.bf.dz;

        % Append to pre-allocated buffer (grow by doubling if full)
        endIdx = nChunkLocs + nDets;
        if endIdx > size(chunkLocs, 1)
            chunkLocs = [chunkLocs; zeros(max(size(chunkLocs,1), nDets), 4)]; %#ok
        end
        frameVec = (globalFrameIdx + iFrame) * ones(nDets, 1);
        chunkLocs(nChunkLocs+1:endIdx, :) = ...
            [x_pos, z_pos, double(maxVals(:)), frameVec];
        nChunkLocs = endIdx;
        t_centroid = t_centroid + toc(tCen);
    end
    chunkLocs = chunkLocs(1:nChunkLocs, :);   % trim pre-allocated tail
    clear IQ_filtered;
    t_det = t_det + toc(tDet);
    
    % O(1) append: just stash the chunk's locs; concatenate once after the loop
    locsPerChunk{iChunk} = chunkLocs;
    allTimestamps(tsFillIdx+1 : tsFillIdx+nFrames) = timestamps;
    tsFillIdx = tsFillIdx + nFrames;
    globalFrameIdx = globalFrameIdx + nFrames;
end

% Concatenate all chunks' localizations in one shot
allLocalizations = vertcat(locsPerChunk{:});
if tsFillIdx < numel(allTimestamps)
    allTimestamps = allTimestamps(1:tsFillIdx);
end

% --- Track within this block ---
tTrack = tic;
if ~isempty(allLocalizations)
    blockTracks = track_microbubbles(allLocalizations, config.track, allTimestamps);
else
    blockTracks = {};
end
t_track = toc(tTrack);

% --- Package result ---
blockResult.localizations  = allLocalizations;
blockResult.tracks         = blockTracks;
blockResult.timestamps     = allTimestamps;
blockResult.numFrames      = globalFrameIdx;
blockResult.numLocs        = size(allLocalizations, 1);
blockResult.numTracks      = numel(blockTracks);
blockResult.processingTime = toc(tBlock);

fprintf('    -> %d locs, %d tracks (%.1f sec)\n', ...
    blockResult.numLocs, blockResult.numTracks, blockResult.processingTime);

% Phase-by-phase timing breakdown
total = blockResult.processingTime;
other = total - (t_load + t_bf + t_mc + t_bg + t_svd + t_det + t_track);
fprintf('    Phase timing: load=%.1fs  bf=%.1fs  mc=%.1fs  bg=%.1fs  svd=%.1fs  det=%.1fs  track=%.1fs  other=%.1fs\n', ...
    t_load, t_bf, t_mc, t_bg, t_svd, t_det, t_track, other);
fprintf('    Det breakdown: detect_microbubbles=%.1fs  centroid=%.1fs\n', ...
    t_detectFn, t_centroid);

blockResult.phaseTiming = struct('load', t_load, 'beamform', t_bf, ...
    'motionCorr', t_mc, 'background', t_bg, 'svd', t_svd, ...
    'detection', t_det, 'tracking', t_track, 'other', other);

end
