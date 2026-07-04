function profile = build_profile(name, overrides)
% BUILD_PROFILE  Named processing profile for the unified IQ engine.
%   profile = build_profile('apr17'|'core', overrides)
%
% PROFILE_APR17 = extracted-verbatim apr17 reproduction (every fix default-OFF);
% process_IQ_block(build_profile('apr17'), IQ, meta) reproduces the defended
% BATCH_FINAL_ANALYSIS inline path. PROFILE_CORE flips the headline defaults.
%
% Runtime-pinned fields (det/track exact values, svd cutoff, lambda, grid) are
% marked "R" in the spec divergence table: the literals here match batch_config.mat,
% but a driver may overwrite them from the loaded cfg/meta for an exact gate.
%
% Headline reanalysis example:
%   build_profile('apr17', struct('localize',struct('mode','centroid'), ...
%       'combineMode','merge', 'svd',struct('seed',12345), ...
%       'det',struct('method','fixed')))
%
% See PHASE2A_BUILD_SPEC_v2.md sec 3.3.

    if nargin < 2, overrides = struct(); end

    switch lower(name)
    case 'apr17'
        % ---- beamform half (BATCH) ----
        profile.useCUDA            = true;          % BATCH:88
        profile.useGPU             = true;          % BATCH:89
        profile.zeroOnly           = false;         % all-5-angle compound
        profile.blankSteering      = false;
        profile.motionCorrection.enable = false;
        profile.bgMeanIQ           = [];            % apr17 does NO bg-IQ subtraction
        profile.sosOverride        = [];            % meta.c = XML ~1480 (BATCH:192)
        profile.reloadMetaPerBlock = true;          % MANDATORY (BATCH:1541-1547)
        profile.chunkSize          = inf;           % whole-block (BATCH:1549)
        profile.timestampSource    = 'frameRate';   % synthesized (BATCH:585,608)
        % ---- SVD ----
        profile.svd.cutoffLow      = 8;             % svdDefault (BATCH:72); LAT-ULM sweep [5,8,12]
        profile.svd.cutoffHigh     = [];            % BATCH:563
        profile.svd.seed           = [];            % UNSEEDED = as-published
        profile.useRSVD            = true;          % BATCH:563
        % ---- detect ----
        profile.det.method         = 'threshold';   % edge-adaptive (BATCH:77)
        profile.det.threshold      = 5;             % BATCH:78
        profile.det.roiSize_px     = 5;             % BATCH:79 (moot for pixel localize)
        profile.det.fixedThresh    = [];            % corrected detector OFF
        profile.det.minSep_mm      = 0.069;         % round(lambda_mm,3) (BATCH:214); pin from meta
        % ---- localize ----
        profile.localize.mode      = 'pixel';       % raw-pixel (BATCH:574-576)
        % ---- track ----
        profile.track.maxDisp_mm        = 0.070;    % BATCH:80
        profile.track.maxGapFrames      = 5;        % BATCH:81
        profile.track.minTrackLength    = 5;        % BATCH:82
        profile.track.kalman.processNoise = 0.010;  % BATCH:83
        profile.track.kalman.measNoise    = 0.080;  % BATCH:84
        profile.track.gapByFrameDelta   = false;    % +1 aging (reproduce retrack)
        profile.track.minMeanAmp        = [];       % QC OFF
        profile.track.minStraightness   = [];       % QC OFF
        profile.track.maxGapInTrack     = [];       % QC OFF
        % ---- combine (driver) ----
        profile.combineMode        = 'retrack';     % BATCH:605-611
        % ---- count axis (driver) ----
        profile.count.metric          = 'tracks';
        profile.count.combineMode     = 'retrack';
        profile.count.allBlocks       = false;      % single = b2 only
        profile.count.bgFloorSubtract = false;
        profile.count.bgRate          = [];
        profile.count.normalize       = 'none';
        % ---- QUASAR ----
        profile.stages.quasar      = true;
        profile.quasar.lambda      = 0.10;          % lambdaDefault (BATCH:71)
        profile.quasar.grid        = 'native';
        profile.quasar.normalize   = false;
        profile.quasar.useNNLS     = false;         % unconstrained-CG-then-clip
        % ---- Doppler ----
        profile.stages.doppler     = false;
        profile.doppler.method     = 'kasai';       % SIGNED lag-1 (BATCH:674-675)
        % ---- grid ----
        profile.gridSpacing        = 'lambda/2';    % BATCH:200-201

    case 'core'
        profile = build_profile('apr17');           % start from apr17, flip headline defaults
        profile.sosOverride        = 1540;          % LAT_ULM:84
        profile.zeroOnly           = true;          % LAT_ULM:59
        profile.blankSteering      = true;          % LAT_ULM:66
        profile.timestampSource    = 'hardware';    % CORE :111
        profile.useRSVD            = false;         % exact SVD (LAT_ULM:48)
        profile.svd.cutoffLow      = 5;             % LAT_ULM:97
        profile.svd.seed           = 12345;         % headline determinism
        profile.det.roiSize_px     = 7;             % LAT_ULM:104
        profile.localize.mode      = 'centroid';    % headline
        profile.track.maxDisp_mm   = 0.500;         % LAT_ULM:107
        profile.track.maxGapFrames = 3;             % LAT_ULM:108
        profile.track.gapByFrameDelta = true;       % headline fix
        profile.combineMode        = 'merge';       % headline fix
        profile.count.combineMode  = 'merge';
        profile.count.allBlocks    = true;
        profile.count.bgFloorSubtract = true;
        profile.quasar.grid        = 'sr';          % super-res grid
        profile.quasar.normalize   = true;
        profile.doppler.method     = 'decorrelation';
        profile.gridSpacing        = 'lambda/5,lambda/10';  % LAT_ULM:609-610

    otherwise
        error('build_profile:unknownName', 'Unknown profile "%s" (expected apr17|core).', name);
    end

    % ---- recursive overlay of user overrides (supports nested structs) ----
    profile = local_apply_overrides(profile, overrides);
end

function s = local_apply_overrides(s, o)
    if isempty(o) || ~isstruct(o), return; end
    f = fieldnames(o);
    for i = 1:numel(f)
        if isstruct(o.(f{i})) && isfield(s, f{i}) && isstruct(s.(f{i}))
            s.(f{i}) = local_apply_overrides(s.(f{i}), o.(f{i}));
        else
            s.(f{i}) = o.(f{i});
        end
    end
end
