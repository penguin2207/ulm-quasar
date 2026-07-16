function detections = detect_microbubbles(frame, params, dx_mm, dz_mm)
% DETECT_MICROBUBBLES  Detect microbubble candidates in a single frame.
%
% Implements two detection approaches:
%   'threshold' - Adaptive thresholding (N * noise_std) with non-maximum
%                 suppression. Simple and fast, suitable when PSF is
%                 well-characterized.
%   'ncc'       - Normalized cross-correlation with a model PSF, followed
%                 by peak detection. More robust to varying signal levels.
%
% Reference: Ackermann D, Schmitz G. "Detection and tracking of multiple
%            microbubbles in ultrasound B-mode images."
%            IEEE Trans UFFC. 2016;63(1):72-82.
%
% Inputs:
%   frame    - [nZ x nX] envelope image (magnitude of IQ)
%   params   - Struct with fields:
%              .method      - 'threshold' or 'ncc'
%              .threshold   - Detection threshold (noise std multiplier)
%              .minSep_mm   - Minimum separation between detections [mm]
%              .roiSize_px  - ROI size for localization [px]
%   dx_mm    - Lateral pixel size [mm]
%   dz_mm    - Axial pixel size [mm]
%
% Output:
%   detections - [N x 2] array of [row, col] pixel indices of detections

[nZ, nX] = size(frame);

switch lower(params.method)
    case 'threshold'
        % Estimate noise from frame edges (regions unlikely to contain bubbles)
        edgePix = min(10, floor(min(nZ, nX) / 4));
        noiseRegion = [frame(1:edgePix, :); frame(end-edgePix+1:end, :)];
        noiseSides  = [frame(edgePix+1:end-edgePix, 1:edgePix), ...
                       frame(edgePix+1:end-edgePix, end-edgePix+1:end)];
        noiseAll    = [noiseRegion(:); noiseSides(:)];
        noiseStd  = std(noiseAll);
        noiseMean = mean(noiseAll);
        
        % Threshold
        thresh = noiseMean + params.threshold * noiseStd;
        
        % Find local maxima above threshold
        detections = find_local_maxima(frame, thresh, params.minSep_mm, dx_mm, dz_mm);
        
    case 'ncc'
        % Build model PSF (2D Gaussian)
        axialSigma_px = (params.psf.axial_mm / dz_mm) / (2 * sqrt(2 * log(2)));
        latSigma_px   = (params.psf.lateral_mm / dx_mm) / (2 * sqrt(2 * log(2)));
        
        halfSize = ceil(3 * max(axialSigma_px, latSigma_px));
        [xx, zz] = meshgrid(-halfSize:halfSize, -halfSize:halfSize);
        psf = exp(-0.5 * (xx.^2 / latSigma_px^2 + zz.^2 / axialSigma_px^2));
        psf = psf / norm(psf(:));
        
        % Normalized cross-correlation
        nccMap = normxcorr2(psf, frame);
        
        % Trim padding from normxcorr2
        pad = halfSize;
        nccMap = nccMap(pad+1:end-pad, pad+1:end-pad);
        
        % Detect peaks
        noiseStd = std(nccMap(:));
        thresh = params.threshold * noiseStd;
        detections = find_local_maxima(nccMap, thresh, params.minSep_mm, dx_mm, dz_mm);

    case 'fixed'
        % CORRECTED DETECTOR. Fixed absolute envelope threshold in envelope-amplitude
        % units, so the floor does NOT move with concentration. This is the fix for the
        % 'threshold' method's real defect (edge-noise * N is signal-dependent: quiet
        % low-conc frames lower the threshold and admit clutter, making the raw count
        % anti-correlated with signal at the low end). See FINDINGS_2026_06_11_lowconc_floor.md.
        %
        % THE THRESHOLD DOES *NOT* YIELD ~0 DETECTIONS ON BUBBLE-FREE Bg BLOCKS. An earlier
        % version of this comment claimed it did, citing calibrate_bg_floor.m; both are wrong
        % and calibrate_bg_floor.m is dead code the runner explicitly rejects. Measured Bg
        % false-alarm rate at the shipped thresholds is 12-18 loc/frame in all six
        % (dataset x domain) combinations, ~1000x the nominal 0.02 tolerance.
        %
        % THE COUNT IS NONETHELESS UNBIASED, because the Bg PEDESTAL SUBTRACTION downstream
        % removes those false alarms exactly: the noise floor is concentration-INDEPENDENT
        % (measured beta = -0.016 vs conc, rung/BGTF 1.04-1.10x across a 147x range), so the
        % bubble-free rate is an unbiased estimate of the rung blocks' false-alarm rate.
        % The threshold sitting inside the noise costs PRECISION, not ACCURACY.
        %
        % DO NOT "FIX" THIS BY RAISING THE THRESHOLD to a false-alarm-free point. Doing so
        % counts only the bright subset, and the bright fraction GROWS with concentration
        % (8.3% at the lowest rung -> 51.1% at the highest, via constructive interference
        % between overlapping bubbles), so a high threshold has RELAXING selectivity as
        % concentration rises. That inflates the top of the ladder and fabricates a slope
        % increase: measured 0.474 -> 0.879, entirely accounted for by the bright-fraction
        % slope of log10(51.1/8.3)/log10(147) = +0.364. A low threshold plus an unbiased
        % pedestal subtraction is CORRECT; a high threshold with no subtraction is
        % amplitude-selected. The threshold is a bias/variance choice, not an error.
        thresh = params.fixedThresh;
        detections = find_local_maxima(frame, thresh, params.minSep_mm, dx_mm, dz_mm);

    otherwise
        error('Unknown detection method: %s', params.method);
end

end


function peaks = find_local_maxima(img, threshold, minSep_mm, dx_mm, dz_mm)
% FIND_LOCAL_MAXIMA  Find local maxima with non-maximum suppression.
%
% Returns pixel coordinates of local maxima above threshold, enforcing
% minimum separation between detections.

[nZ, nX] = size(img);

% Minimum separation in pixels
minSepZ = ceil(minSep_mm / dz_mm);
minSepX = ceil(minSep_mm / dx_mm);

% Find all pixels above threshold
candidates = img > threshold;

% 3x3 local maximum test via separable moving-max (movmax is base MATLAB
% and much faster than the 8-comparison loop: one pass along each axis
% instead of eight shifted-image allocations and ANDs).
%
% movmax default 'shrink' endpoint handling is equivalent to -inf padding
% for a max filter: out-of-range samples simply don't affect the max.
max3x3 = movmax(movmax(img, 3, 1), 3, 2);
localMax = (img == max3x3);

candidates = candidates & localMax;

% Extract candidate list sorted by intensity (brightest first)
[rows, cols] = find(candidates);
if isempty(rows)
    peaks = zeros(0, 2);
    return;
end
vals = img(sub2ind(size(img), rows, cols));
[~, sortIdx] = sort(vals, 'descend');
rows = rows(sortIdx);
cols = cols(sortIdx);

% Non-maximum suppression: grid-based O(N) approach
% Assign each candidate to a grid cell; keep only the brightest per cell
gridZ = ceil(rows / max(minSepZ, 1));
gridX = ceil(cols / max(minSepX, 1));
nGZ = max(gridZ); nGX = max(gridX);
occupied = false(nGZ, nGX);
keep = false(numel(rows), 1);
for i = 1:numel(rows)
    gz = gridZ(i); gx = gridX(i);
    % Check this cell and neighbors
    g1z = max(1, gz-1); g2z = min(nGZ, gz+1);
    g1x = max(1, gx-1); g2x = min(nGX, gx+1);
    if ~any(occupied(g1z:g2z, g1x:g2x), 'all')
        keep(i) = true;
        occupied(gz, gx) = true;
    end
end

peaks = [rows(keep), cols(keep)];

end
