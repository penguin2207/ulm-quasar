function qc = assess_block_quality(dataFolder, baseFilename, modeName, ...
    config, anglePairs, elemPos_mm, delayTables, xGrid, zGrid, nX, nZ, ...
    fs_MHz, c, depthOffset_mm, bgMeanIQ)
% ASSESS_BLOCK_QUALITY  Quick screen of a VADA block for bubble presence.
%
% Beamforms a small sample of frames (~50), computes temporal std and
% detection metrics, and returns a quality struct for block discrimination.
%
% Metrics computed:
%   stdContrast   - Peak/median of temporal std (high = localized motion)
%   locsStrict    - Detections/frame at strict threshold (>0 = likely bubbles)
%   locsMild      - Detections/frame at mild threshold (comparison)
%   snrDrop       - SNR difference between SVD cut=2 and cut=20
%                   (>2 dB means signal in mid-range SVs = bubbles)
%   hasBubbles    - Boolean: recommended include/exclude

tStart = tic;

qc = struct();
qc.filename = baseFilename;
qc.stdContrast = 0;
qc.locsStrict = 0;
qc.locsMild = 0;
qc.snrDrop = 0;
qc.hasBubbles = false;
qc.assessTime = 0;

% Sample size: 50 compound frames is enough for quick assessment
numSampleFrames = 50;
eventsPerFrame = config.eventsPerFrame;
numSampleEvents = numSampleFrames * eventsPerFrame;

try
    % Check total events available
    [VadaTest, ~, ~, BlockConfig] = VsiVadaDataRead(dataFolder, baseFilename, ...
        1:eventsPerFrame, modeName);
    numTotalEvents = numel(BlockConfig.PulseSequences(1).Events);
    clear VadaTest;
    
    numSampleEvents = min(numSampleEvents, numTotalEvents);
    numSampleFrames = floor(numSampleEvents / eventsPerFrame);
    
    if numSampleFrames < 10
        qc.assessTime = toc(tStart);
        return;  % Too few frames to assess
    end
    
    % Load sample events
    [VadaSample, ~, ~, ~] = VsiVadaDataRead(dataFolder, baseFilename, ...
        1:numSampleEvents, modeName);
    
    % Beamform to complex IQ
    iqStack = zeros(nZ, nX, numSampleFrames, 'single');
    
    if config.zeroOnly
        for iFrame = 1:numSampleFrames
            evIdx = (iFrame-1) * eventsPerFrame + config.zeroEventIdx;
            rfData = single(VadaSample(evIdx).Data);
            rxPos_mm = elemPos_mm(anglePairs(config.zeroAngleIdx).rxElements);
            bfImg = beamform_planewave_gpu(rfData, rxPos_mm, ...
                anglePairs(config.zeroAngleIdx).angle, ...
                xGrid, zGrid, fs_MHz, c, depthOffset_mm, delayTables{1});
            iqStack(:,:,iFrame) = single(bfImg);
        end
    else
        for iFrame = 1:numSampleFrames
            baseEvt = (iFrame-1) * eventsPerFrame;
            compImg = complex(zeros(nZ, nX, 'single'));
            for a = 1:config.numAngles
                if config.hasPI
                    rfData = single(VadaSample(baseEvt + anglePairs(a).posIdx).Data) + ...
                             single(VadaSample(baseEvt + anglePairs(a).negIdx).Data);
                else
                    rfData = single(VadaSample(baseEvt + anglePairs(a).posIdx).Data);
                end
                rxPos_mm = elemPos_mm(anglePairs(a).rxElements);
                bfImg = beamform_planewave_gpu(rfData, rxPos_mm, anglePairs(a).angle, ...
                    xGrid, zGrid, fs_MHz, c, depthOffset_mm, delayTables{a});
                compImg = compImg + single(bfImg);
            end
            iqStack(:,:,iFrame) = compImg;
        end
    end
    clear VadaSample;
    
    % Background subtraction
    if ~isempty(bgMeanIQ)
        for iFrame = 1:numSampleFrames
            iqStack(:,:,iFrame) = iqStack(:,:,iFrame) - bgMeanIQ;
        end
    end
    
    % Metric 1: Temporal std contrast (peak / median)
    complexStd = abs(std(iqStack, 0, 3));
    qc.stdContrast = max(complexStd(:)) / (median(complexStd(:)) + eps);
    
    % SVD
    Casorati = reshape(iqStack, nZ*nX, numSampleFrames);
    [U, S, V] = svd(single(Casorati), 'econ');
    singVals = diag(S);
    clear Casorati;
    
    % Metric 2: SNR drop between low and high SVD cutoff
    % (real bubbles live in mid-range SVs, so removing them drops SNR)
    cut_lo = min(2, numSampleFrames-1);
    cut_hi = min(20, numSampleFrames-1);
    
    filt_lo = U(:, cut_lo+1:end) * S(cut_lo+1:end, cut_lo+1:end) * V(:, cut_lo+1:end)';
    env_lo = mean(abs(reshape(filt_lo, nZ, nX, numSampleFrames)), 3);
    noise_lo = env_lo(1:min(20,nZ), :);
    snr_lo = 20*log10(max(env_lo(:)) / (std(noise_lo(:)) + eps));
    
    filt_hi = U(:, cut_hi+1:end) * S(cut_hi+1:end, cut_hi+1:end) * V(:, cut_hi+1:end)';
    env_hi = mean(abs(reshape(filt_hi, nZ, nX, numSampleFrames)), 3);
    noise_hi = env_hi(1:min(20,nZ), :);
    snr_hi = 20*log10(max(env_hi(:)) / (std(noise_hi(:)) + eps));
    
    qc.snrDrop = snr_lo - snr_hi;
    
    % Metric 3: Detection counts at strict and mild thresholds
    % Use SVD cutoff=2 filtered data
    filtStack = reshape(filt_lo, nZ, nX, numSampleFrames);
    dx = xGrid(2) - xGrid(1);
    dz = zGrid(2) - zGrid(1);
    
    detParams.method = 'threshold';
    detParams.minSep_mm = 0.200;
    detParams.roiSize_px = 7;
    
    % Strict: threshold=8 (only real signal survives)
    detParams.threshold = 8;
    strictLocs = 0;
    for f = 1:numSampleFrames
        dets = detect_microbubbles(abs(filtStack(:,:,f)), detParams, dx, dz);
        strictLocs = strictLocs + size(dets, 1);
    end
    qc.locsStrict = strictLocs / numSampleFrames;
    
    % Mild: threshold=5
    detParams.threshold = 5;
    mildLocs = 0;
    for f = 1:numSampleFrames
        dets = detect_microbubbles(abs(filtStack(:,:,f)), detParams, dx, dz);
        mildLocs = mildLocs + size(dets, 1);
    end
    qc.locsMild = mildLocs / numSampleFrames;
    
    % Decision: block has bubbles if any of these criteria met
    qc.hasBubbles = (qc.locsStrict > 0.05) || ...   % Any strict detections
                    (qc.snrDrop > 2.0) || ...         % Signal in mid-range SVs
                    (qc.stdContrast > 3.0);            % Localized temporal variance
    
catch ME
    fprintf('      QC error: %s\n', ME.message);
end

qc.assessTime = toc(tStart);

end
