function [shifts, diagnostics] = estimate_tissue_motion(IQ_stack, params)
% ESTIMATE_TISSUE_MOTION  Estimate global rigid motion from IQ image stack.
%
% Uses phase correlation on the tissue envelope (magnitude of complex IQ)
% to estimate frame-to-frame axial and lateral shifts caused by cardiac
% or respiratory motion. Sub-pixel precision via parabolic peak fitting.
%
% This should be called BEFORE SVD clutter filtering, since tissue motion
% degrades SVD separation of tissue and blood signal components.
%
% Reference: Demene et al., IEEE TMI 2015;34(11):2271-2285
%
% Inputs:
%   IQ_stack  - [nZ x nX x nFrames] complex IQ image stack
%   params    - Struct with fields:
%       .method    - 'phase_corr' (default)
%       .refType   - 'rolling' (default) or 'first'
%       .refWindow - Rolling window half-width in frames (default: 10)
%       .maxShift  - Max expected shift in pixels (default: 5)
%
% Outputs:
%   shifts      - [nFrames x 2] estimated shifts [dz, dx] in pixels
%                  (fractional for sub-pixel). Frame 1 is always [0, 0].
%   diagnostics - Struct with:
%       .maxDisp_px   - Maximum displacement magnitude [pixels]
%       .meanDisp_px  - Mean displacement magnitude [pixels]
%       .peakCorr     - [nFrames x 1] peak correlation values (quality)

[nZ, nX, nFrames] = size(IQ_stack);

% Defaults
if nargin < 2, params = struct(); end
if ~isfield(params, 'method'),    params.method    = 'phase_corr'; end
if ~isfield(params, 'refType'),   params.refType   = 'rolling';    end
if ~isfield(params, 'refWindow'), params.refWindow = 10;           end
if ~isfield(params, 'maxShift'),  params.maxShift  = 5;            end

shifts = zeros(nFrames, 2);   % [dz, dx]
peakCorr = ones(nFrames, 1);  % Correlation quality metric

% Compute tissue envelope (magnitude)
envelope = abs(IQ_stack);

% Precompute reference for 'first' mode
if strcmpi(params.refType, 'first')
    refImg = envelope(:,:,1);
end

maxS = params.maxShift;

for iFrame = 2:nFrames
    % Build reference image
    if strcmpi(params.refType, 'rolling')
        % Rolling mean centered on current frame
        w = params.refWindow;
        rStart = max(1, iFrame - w);
        rEnd   = min(nFrames, iFrame + w);
        % Exclude current frame from reference
        refIdx = [rStart:iFrame-1, iFrame+1:rEnd];
        if isempty(refIdx), refIdx = max(1, iFrame-1); end
        refImg = mean(envelope(:,:,refIdx), 3);
    end
    % else 'first' mode: refImg already set

    curImg = envelope(:,:,iFrame);

    % Phase correlation
    F_ref = fft2(refImg);
    F_cur = fft2(curImg);

    % Cross-power spectrum
    crossPower = F_cur .* conj(F_ref);
    crossPower = crossPower ./ (abs(crossPower) + eps);

    % Inverse FFT gives phase correlation surface
    corrSurf = real(ifft2(crossPower));

    % Shift so that zero-lag is at center
    corrSurf = fftshift(corrSurf);
    centerZ = ceil((nZ + 1) / 2);
    centerX = ceil((nX + 1) / 2);

    % Restrict search to maxShift region around center
    zLo = max(1,  centerZ - maxS);
    zHi = min(nZ, centerZ + maxS);
    xLo = max(1,  centerX - maxS);
    xHi = min(nX, centerX + maxS);

    searchRegion = corrSurf(zLo:zHi, xLo:xHi);

    % Find integer peak
    [maxVal, linIdx] = max(searchRegion(:));
    [pz, px] = ind2sub(size(searchRegion), linIdx);

    % Convert to global coordinates
    peakZ = pz + zLo - 1;
    peakX = px + xLo - 1;

    % Sub-pixel refinement via parabolic fit (3-point)
    dz_sub = 0;
    dx_sub = 0;

    if peakZ > 1 && peakZ < nZ
        vPrev = corrSurf(peakZ-1, peakX);
        vCurr = corrSurf(peakZ,   peakX);
        vNext = corrSurf(peakZ+1, peakX);
        denom = 2 * (2*vCurr - vPrev - vNext);
        if abs(denom) > eps
            dz_sub = (vPrev - vNext) / denom;
        end
    end

    if peakX > 1 && peakX < nX
        vPrev = corrSurf(peakZ, peakX-1);
        vCurr = corrSurf(peakZ, peakX);
        vNext = corrSurf(peakZ, peakX+1);
        denom = 2 * (2*vCurr - vPrev - vNext);
        if abs(denom) > eps
            dx_sub = (vPrev - vNext) / denom;
        end
    end

    % Total shift relative to reference (negative = moved in positive direction)
    shifts(iFrame, 1) = (peakZ - centerZ) + dz_sub;  % dz
    shifts(iFrame, 2) = (peakX - centerX) + dx_sub;  % dx
    peakCorr(iFrame) = maxVal;
end

% Diagnostics
dispMag = sqrt(shifts(:,1).^2 + shifts(:,2).^2);
diagnostics.maxDisp_px  = max(dispMag);
diagnostics.meanDisp_px = mean(dispMag);
diagnostics.peakCorr    = peakCorr;
diagnostics.shifts      = shifts;

end
