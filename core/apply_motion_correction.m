function IQ_corrected = apply_motion_correction(IQ_stack, shifts, useGPU)
% APPLY_MOTION_CORRECTION  Sub-pixel rigid motion correction via Fourier shift.
%
% Applies the inverse of estimated tissue motion to realign complex IQ
% frames. Uses the Fourier shift theorem for sub-pixel precision while
% perfectly preserving phase information (critical for SVD and PI).
%
% Fourier shift theorem:
%   Shifted image = ifft2( fft2(image) .* exp(-j*2*pi*(fz*dz + fx*dx)) )
%
% Must be applied to COMPLEX IQ data (not magnitude), since downstream
% SVD filtering and pulse inversion depend on phase coherence.
%
% Inputs:
%   IQ_stack  - [nZ x nX x nFrames] complex IQ image stack
%   shifts    - [nFrames x 2] shifts to correct [dz, dx] in pixels
%               (as returned by estimate_tissue_motion)
%   useGPU    - Boolean: use GPU acceleration for fft2/ifft2
%
% Output:
%   IQ_corrected - [nZ x nX x nFrames] motion-corrected complex IQ stack

[nZ, nX, nFrames] = size(IQ_stack);

if nargin < 3, useGPU = false; end

IQ_corrected = zeros(nZ, nX, nFrames, 'like', IQ_stack);

% Build frequency grids (normalized: 0 to 1 in each dimension)
[fx, fz] = meshgrid( ...
    ([0:nX-1] - floor(nX/2)) / nX, ...
    ([0:nZ-1] - floor(nZ/2)) / nZ);

% ifftshift to align with fft2 output ordering
fz = ifftshift(fz);
fx = ifftshift(fx);

if useGPU
    fz = gpuArray(single(fz));
    fx = gpuArray(single(fx));
end

for iFrame = 1:nFrames
    dz = shifts(iFrame, 1);
    dx = shifts(iFrame, 2);

    if abs(dz) < 1e-6 && abs(dx) < 1e-6
        % No shift needed
        IQ_corrected(:,:,iFrame) = IQ_stack(:,:,iFrame);
        continue;
    end

    frame = IQ_stack(:,:,iFrame);
    if useGPU
        frame = gpuArray(single(frame));
    end

    % Fourier shift: apply linear phase ramp in frequency domain
    % Negative shift to UNDO the estimated motion
    phaseShift = exp(-1j * 2 * pi * (fz * (-dz) + fx * (-dx)));

    F_frame = fft2(frame);
    F_shifted = F_frame .* phaseShift;
    corrected = ifft2(F_shifted);

    if useGPU
        corrected = gather(corrected);
    end

    IQ_corrected(:,:,iFrame) = corrected;
end

end
