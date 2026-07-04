function [velMap, decorMap] = estimate_decorrelation_velocity(IQ_ensemble, ...
    frameRate, lambda_mm, srFactor, useGPU)
% ESTIMATE_DECORRELATION_VELOCITY  Velocity from speckle decorrelation.
%
% Within a SUSHI ensemble, flowing microbubbles cause the complex IQ
% signal to decorrelate over time. The decorrelation rate is proportional
% to velocity: faster flow = faster decorrelation.
%
% Method:
%   1. Compute normalized lag-1 autocorrelation of complex IQ per pixel:
%      R(1) = mean(IQ(t) .* conj(IQ(t+1))) / mean(|IQ(t)|^2)
%   2. Decorrelation coefficient: D = 1 - |R(1)|
%   3. Velocity estimate: v = lambda / (2*pi*tau)
%      where tau = -1 / (frameRate * log(|R(1)|))
%
% Inputs:
%   IQ_ensemble  - [nZ x nX x nFrames] complex IQ data for one ensemble
%   frameRate    - Compound frame rate [Hz]
%   lambda_mm    - Wavelength [mm]
%   srFactor     - Super-resolution upsampling factor
%   useGPU       - Use GPU (default false)
%
% Outputs:
%   velMap    - [nZ*srFactor x nX*srFactor] velocity map [mm/s]
%   decorMap  - [nZ*srFactor x nX*srFactor] decorrelation coefficient (0-1)

if nargin < 5, useGPU = false; end

[nZ, nX, nFrames] = size(IQ_ensemble);

if nFrames < 3
    velMap = zeros(nZ * srFactor, nX * srFactor, 'single');
    decorMap = velMap;
    return;
end

% Compute lag-1 complex autocorrelation per pixel
% R1 = sum(IQ(t) .* conj(IQ(t+1))) / sum(|IQ(t)|^2)
if useGPU
    IQ = gpuArray(single(IQ_ensemble));
else
    IQ = single(IQ_ensemble);
end

% Numerator: cross-correlation at lag 1
R1_num = sum(IQ(:,:,1:end-1) .* conj(IQ(:,:,2:end)), 3);

% Denominator: power (use geometric mean of both lags for stability)
P0 = sum(abs(IQ(:,:,1:end-1)).^2, 3);
P1 = sum(abs(IQ(:,:,2:end)).^2, 3);
R1_den = sqrt(P0 .* P1) + eps;

% Normalized lag-1 autocorrelation
R1 = R1_num ./ R1_den;

if useGPU
    R1 = gather(R1);
end

% Decorrelation coefficient: 0 = static, 1 = fully decorrelated
absR1 = abs(R1);
decorCoeff = 1 - absR1;

% Velocity from decorrelation time
% |R(1)| = exp(-dt/tau), so tau = -dt / log(|R(1)|)
% v = lambda / (2*pi*tau)
dt = 1 / frameRate;
logR1 = log(max(absR1, 1e-6));  % Clamp to avoid log(0)
tau = -dt ./ logR1;
velocity = lambda_mm ./ (2 * pi * tau);

% Clamp unreasonable velocities (static tissue gives |R1|~1, tau~inf, v~0)
velocity(absR1 > 0.999) = 0;    % Nearly static
velocity(absR1 < 0.01) = NaN;   % Pure noise (fully decorrelated)
velocity = max(velocity, 0);     % No negative velocities

% Upsample to SR grid
velMap = single(imresize(velocity, [nZ * srFactor, nX * srFactor], 'bicubic'));
decorMap = single(imresize(decorCoeff, [nZ * srFactor, nX * srFactor], 'bicubic'));

% Clean up NaN propagation from upsampling
velMap(isnan(velMap)) = 0;
decorMap(isnan(decorMap)) = 0;

end
