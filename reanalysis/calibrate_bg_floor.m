function cal = calibrate_bg_floor(bgIQ, profile, meta, opts)
% CALIBRATE_BG_FLOOR  Fixed detection threshold + Bg track-floor from bubble-free blocks.
%
%   cal = calibrate_bg_floor(bgIQ, profile, meta, opts)
%
% Calibrates profile.det.fixedThresh (for det.method='fixed') so the bubble-free Bg
% blocks yield ~0 detections. Per PHASE2A_BUILD_SPEC_v2.md Part B / Risk R-DET-F1:
% the threshold is swept in LOC-RATE space (running the full detector incl 3x3 NMS +
% minSep), NOT chosen by pixel-percentile (NMS/minSep make a pixel-percentile
% non-monotonic with loc/frame). The Bg blocks are the calibration instrument; see
% FINDINGS_2026_06_11_lowconc_floor.md (Bg ~5 false loc/frame on the OLD adaptive detector).
%
% bgIQ    : cell {IQ_bg1, IQ_bg2, ...} of loaded bubble-free IQ stacks.
% profile : build_profile output; uses .svd (cutoff+seed), .det, .track, .useGPU.
% meta    : carries .dx .dz .xGrid .zGrid .nZ .nX .frameRate_Hz.
% opts    : .tolLocPerFrame (default 0.05), .nThr (default 25), .nFrCal (default 600).
%
% Returns:
%   cal.fixedThresh  - chosen absolute envelope threshold (inject into profile.det.fixedThresh)
%   cal.bgLocRate    - [1 x nBg] loc/frame at the chosen threshold (<= tol if achievable)
%   cal.bgTrackRate  - [1 x nBg] QC-track/frame at the chosen threshold (the floor to subtract)
%   cal.sweep        - struct(thrs, locRate[nThr x nBg], tol, cutoff, chosenIdx)

    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'tolLocPerFrame'), opts.tolLocPerFrame = 8; end   % SPARSE loc/frame cap: keep tracking meaningful (not a clutter storm)
    if ~isfield(opts,'nThr'),    opts.nThr    = 25;  end
    if ~isfield(opts,'nFrCal'),  opts.nFrCal  = 600; end

    useGPU = profile.useGPU; seed = profile.svd.seed; cut = profile.svd.cutoffLow;
    nBg = numel(bgIQ);

    % --- SVD-filter each Bg block (seeded, same cutoff as the conc rungs) ---
    filt = cell(1, nBg); nFr = zeros(1, nBg); envPool = [];
    for b = 1:nBg
        IQ = bgIQ{b};
        if size(IQ,3) > opts.nFrCal, IQ = IQ(:,:,1:opts.nFrCal); end
        filt{b} = svd_clutter_filter_rsvd(IQ, cut, [], useGPU, seed);
        nFr(b)  = size(filt{b}, 3);
        e = abs(filt{b});
        step = max(1, round(numel(e)/2e4));
        envPool = [envPool; e(1:step:end)']; %#ok<AGROW> subsample for the bracket
    end

    % --- bracket the sweep (percentile only SEEDS the range, not the value) ---
    loThr = prctile(envPool, 50);
    hiThr = prctile(envPool, 99.99);
    thrs  = linspace(loThr, hiThr, opts.nThr);

    if ~isfield(opts,'tolTrackPerFrame'), opts.tolTrackPerFrame = 0.002; end  % ~<=1 track / 500 frames
    if ~isfield(opts,'locTrackCap'),      opts.locTrackCap      = 200;     end  % skip tracking clutter storms

    % Sweep thr. COUNT locs cheaply at every thr; only COLLECT+TRACK+QC where the
    % loc-rate is already near the floor (<= locTrackCap). The COUNT axis counts
    % QC-tracks, so the track-rate is the floor that matters: the min-length+kinematic
    % filter kills clutter over a wide thr range, so we pick the SMALLEST thr that keeps
    % Bg tracks ~0 (max low-conc sensitivity), not the (over-strict) loc-rate=0 point.
    % (Tracking millions of clutter detections at low thr is pointless and O(n^2).)
    det = profile.det; det.method = 'fixed';
    locRate   = zeros(opts.nThr, nBg);
    trackRate = inf(opts.nThr, nBg);
    for it = 1:opts.nThr
        det.fixedThresh = thrs(it);
        for b = 1:nBg
            tot = 0;
            for fr = 1:nFr(b)
                pix = detect_microbubbles(abs(filt{b}(:,:,fr)), det, meta.dx, meta.dz);
                tot = tot + size(pix,1);
            end
            locRate(it,b) = tot / nFr(b);
            if locRate(it,b) == 0
                trackRate(it,b) = 0;
            elseif locRate(it,b) <= opts.locTrackCap
                locs  = local_collect_locs(filt{b}, det, meta);
                ts    = ((0:nFr(b)-1)' * 1000 / meta.frameRate_Hz);
                trk   = track_microbubbles(locs, profile.track, ts);
                trkQC = filter_tracks_quality(trk, profile.track);
                trackRate(it,b) = numel(trkQC) / nFr(b);
            end   % else: loc-rate too high -> trackRate stays inf (above the floor)
        end
    end

    % --- pick smallest thr where ALL Bg blocks are BOTH sparse (loc-rate <= tolLoc) AND
    %     Bg-clean (QC-track-rate <= tolTrack). CONJUNCTIVE: the loc-rate gate is ESSENTIAL.
    %     Track-rate alone lands at a clutter-storm threshold (~55 loc/frame) because the
    %     length filter zeroes Bg tracks long before detections are sparse, so tracking is
    %     no longer a meaningful discriminator. [bug found by verify-reanalysis workflow] ---
    ok  = all(trackRate <= opts.tolTrackPerFrame, 2) & all(locRate <= opts.tolLocPerFrame, 2);
    idx = find(ok, 1, 'first');
    if isempty(idx)
        [~, idx] = min(max(trackRate, [], 2));  % closest if none clears tol
        warning('calibrate_bg_floor:noTol', ...
            'No threshold met both loc(<=%.3g)+track(<=%.3g)/frame tol; using closest track-rate.', ...
            opts.tolLocPerFrame, opts.tolTrackPerFrame);
    end
    cal.fixedThresh = thrs(idx);
    cal.bgLocRate   = locRate(idx,:);
    cal.bgTrackRate = trackRate(idx,:);
    cal.sweep = struct('thrs', thrs, 'locRate', locRate, 'trackRate', trackRate, ...
                       'tolTrack', opts.tolTrackPerFrame, 'cutoff', cut, 'chosenIdx', idx);
end

function locs = local_collect_locs(IQ_filt, det, meta)
    nFr = size(IQ_filt,3); C = cell(nFr,1);
    for fr = 1:nFr
        env = abs(IQ_filt(:,:,fr));
        pix = detect_microbubbles(env, det, meta.dx, meta.dz);
        if isempty(pix), continue; end
        xc = meta.xGrid(pix(:,2)); zc = meta.zGrid(pix(:,1));
        a  = env(sub2ind([meta.nZ meta.nX], pix(:,1), pix(:,2)));
        C{fr} = [xc(:) zc(:) a(:) fr*ones(size(pix,1),1)];
    end
    locs = vertcat(C{:});
end
