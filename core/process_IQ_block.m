function blockResult = process_IQ_block(profile, IQ, meta)
% PROCESS_IQ_BLOCK  Unified IQ->blockResult engine (SVD->detect->localize->track
% + QUASAR/SUSHI + optional Doppler), profile-selected.
%
%   blockResult = process_IQ_block(profile, IQ, meta)
%
% Cut point = process_single_block.m:240 (the clean beamform->IQ boundary).
% Stateless, per-block. PROFILE_APR17 reproduces the BATCH inline path
% (pixel localize, QUASAR on ensemble power); PROFILE_CORE uses sub-pixel centroid.
%
% IQ   : [nZ x nX x nFrames] single complex, PI-sum/angle-compounded CONTRAST domain.
% meta : per-block bundle, MUST carry .dx .dz .xGrid .zGrid .nZ .nX .frameRate_Hz
%        .globalFrameOffset .psf .blockName.
%
% See PHASE2A_BUILD_SPEC_v2.md sec 3.2.

    useGPU = profile.useGPU;
    if ~isfield(meta,'globalFrameOffset') || isempty(meta.globalFrameOffset)
        meta.globalFrameOffset = 0;
    end
    blockResult = struct('filename', meta.blockName);

    % ===== SVD clutter filter (CORE :242-249 / BATCH :563) =====
    % The seed hook fires INSIDE svd_clutter_filter_rsvd before the randn Omega.
    IQ_filt = svd_clutter_filter_rsvd(IQ, profile.svd.cutoffLow, ...
                                      profile.svd.cutoffHigh, useGPU, profile.svd.seed);

    % ===== LAT-ULM stage (always on) =====
    nFr = size(IQ_filt, 3);
    locsCell = cell(nFr, 1);
    for fr = 1:nFr
        env = abs(IQ_filt(:,:,fr));                            % CORE :282 / BATCH :571
        pix = detect_microbubbles(env, profile.det, meta.dx, meta.dz);  % :283 / :572
        if isempty(pix), continue; end

        switch profile.localize.mode
            case 'pixel'                                       % apr17 (BATCH:574-576)
                x_mm = meta.xGrid(pix(:,2)); x_mm = x_mm(:);
                z_mm = meta.zGrid(pix(:,1)); z_mm = z_mm(:);
                amp  = env(sub2ind([meta.nZ meta.nX], pix(:,1), pix(:,2)));  % peak-pixel env
                amp  = amp(:);
            case 'centroid'                                    % headline (psb:288-344)
                [x_mm, z_mm, amp] = local_centroid_localize(env, pix, meta, profile);
            otherwise
                error('process_IQ_block:badLocalize', 'Unknown localize.mode "%s".', profile.localize.mode);
        end
        frameIdx = (meta.globalFrameOffset + fr) * ones(size(pix,1), 1);   % BATCH:577
        locsCell{fr} = [x_mm, z_mm, amp, frameIdx];
    end
    locs = vertcat(locsCell{:});

    ts = ((meta.globalFrameOffset:meta.globalFrameOffset+nFr-1)' * 1000 / meta.frameRate_Hz);  % :585
    % Density guard: an uncalibrated / too-low detector threshold floods locs and
    % makes tracking intractable (track_microbubbles is superlinear in locs, so a
    % 10 h batch would hang). Real ULM is well under a few hundred locs/frame, so
    % skip + flag a block whose density is absurd rather than hang on it.
    maxLPF = 800;
    if isfield(profile,'track') && isfield(profile.track,'maxLocsPerFrame') ...
            && ~isempty(profile.track.maxLocsPerFrame)
        maxLPF = profile.track.maxLocsPerFrame;
    end
    locsPerFrame = size(locs,1) / max(nFr,1);
    blockResult.tooDense = false;
    if locsPerFrame > maxLPF
        warning('process_IQ_block:tooDense', ['%s: %.0f locs/frame exceeds cap %.0f -- ' ...
            'detector threshold too low (uncalibrated?); SKIPPING tracking, counts flagged ' ...
            'for this block.'], meta.blockName, locsPerFrame, maxLPF);
        blockResult.tooDense = true;
        tracks = {};
    elseif ~isempty(locs)
        tracks = track_microbubbles(locs, profile.track, ts);             % :587 / CORE :377
    else
        tracks = {};
    end
    tracks = tracks(:);
    tracksQC = filter_tracks_quality(tracks, profile.track);              % NEW, default pass-through

    blockResult.localizations     = locs;
    blockResult.tracks            = tracks;
    blockResult.numLocs           = size(locs, 1);                        % raw (floor-artifact axis)
    blockResult.numLocsRaw        = size(locs, 1);
    blockResult.numTracks         = numel(tracks);                        % length-filtered
    blockResult.numTracksFiltered = numel(tracks);
    blockResult.numTracksQC       = numel(tracksQC);                      % length + quality
    blockResult.numFrames         = nFr;
    blockResult.timestamps        = ts;
    blockResult.localizeMode      = profile.localize.mode;
    if strcmp(profile.localize.mode, 'centroid')
        blockResult.ampDefinition = 'roiMax';
    else
        blockResult.ampDefinition = 'peakPixelEnv';
    end

    % ===== QUASAR / SUSHI stage (optional; BATCH :461-483) =====
    if profile.stages.quasar
        powerImg = mean(abs(IQ_filt).^2, 3);                              % :465 ensemble power
        [x_s, sInfo] = sushi_sparse_recovery(powerImg, meta.psf, profile.quasar.lambda, ...
            'fista', struct('maxIter',100, 'tol',1e-4, 'nonNeg',true, 'useGPU',useGPU));  % :471-472
        qopts = struct('maxIterCG',50, 'tolCG',1e-6, 'useGPU',useGPU);    % :476-477
        if isfield(profile.quasar,'useNNLS') && profile.quasar.useNNLS
            qopts.useNNLS = true;                                         % requires modified quasar_refit
        end
        [x_q, qInfo] = quasar_refit(powerImg, meta.psf, x_s, qopts);
        blockResult.quasar.sumSushi  = sum(x_s(:));                       % :480
        blockResult.quasar.sumQuasar = sum(x_q(:));                       % :481
        blockResult.quasar.support   = sInfo.support;                    % :482
        blockResult.quasar.ampRatio  = qInfo.amplitudeRatio;             % :483
        blockResult.quasar.x_sushi   = x_s;
        blockResult.quasar.x_quasar  = x_q;
        blockResult.quasar.lambda    = profile.quasar.lambda;
        blockResult.quasar.svdCutoff = profile.svd.cutoffLow;
    end

    % ===== Doppler stage (optional; BATCH :674-675) =====
    if isfield(profile.stages,'doppler') && profile.stages.doppler
        R1 = mean(IQ_filt(:,:,2:end) .* conj(IQ_filt(:,:,1:end-1)), 3);   % Kasai lag-1
        blockResult.doppler.phase  = angle(R1);                          % SIGNED
        blockResult.doppler.method = profile.doppler.method;
    end

    blockResult.svdCutoff = profile.svd.cutoffLow;
end

% ------------------------------------------------------------------------
function [x_mm, z_mm, amp] = local_centroid_localize(frame, dets, meta, profile)
% Vectorized sub-pixel intensity-weighted centroid, faithful to
% process_single_block.m:288-344. amp = ROI max (NOT peak-pixel env).
    halfROI = floor(profile.det.roiSize_px / 2);
    roiSize = 2*halfROI + 1;
    roiPx   = roiSize * roiSize;
    [colGrid, rowGrid] = meshgrid(1:roiSize, 1:roiSize);
    rowVec = single(rowGrid(:));
    colVec = single(colGrid(:));
    nDets  = size(dets, 1);

    frame_padded = padarray(single(frame), [halfROI halfROI], 0, 'both');
    nZp = size(frame_padded, 1);
    det_r_pad = int32(dets(:,1)) + int32(halfROI);
    det_c_pad = int32(dets(:,2)) + int32(halfROI);
    offsets_col = int32(-halfROI:halfROI).';
    offsets_row = reshape(offsets_col, 1, roiSize);
    r_idx = reshape(det_r_pad, 1, 1, nDets) + offsets_col;
    c_idx = reshape(det_c_pad, 1, 1, nDets) + offsets_row;
    lin_idx = r_idx + (c_idx - 1) * int32(nZp);
    rois = frame_padded(lin_idx);
    rois_flat = reshape(rois, roiPx, nDets);

    roi_mins = min(rois_flat, [], 1);
    rois_bg  = max(rois_flat - roi_mins, 0);
    tots  = sum(rois_bg, 1);
    sum_r = rowVec' * rois_bg;
    sum_c = colVec' * rois_bg;
    safe_tots = max(tots, eps('single'));
    subR = sum_r ./ safe_tots;
    subC = sum_c ./ safe_tots;
    zero_mask = (tots == 0);
    center_val = single((roiSize + 1) / 2);
    subR(zero_mask) = center_val;
    subC(zero_mask) = center_val;
    maxVals = max(rois_flat, [], 1);

    x_center = reshape(meta.xGrid(dets(:,2)), nDets, 1);
    z_center = reshape(meta.zGrid(dets(:,1)), nDets, 1);
    x_mm = x_center + (double(subC(:)) - halfROI - 1) * meta.dx;
    z_mm = z_center + (double(subR(:)) - halfROI - 1) * meta.dz;
    amp  = double(maxVals(:));
end
