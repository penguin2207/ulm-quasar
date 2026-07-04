function bfImg = beamform_planewave_gpu(rfData, rxPos_mm, steerAngle_deg, ...
    xGrid, zGrid, fs_MHz, c, depthOffset_mm, delayTable)
% BEAMFORM_PLANEWAVE_GPU  GPU-accelerated delay-and-sum beamforming.
%
% Fully vectorized DAS beamforming using gpuArray. All delay computation
% and interpolation runs on GPU. With an A6000 (48 GB VRAM), entire delay
% tables and RF data fit comfortably in GPU memory.
%
% Reference: Montaldo et al., IEEE TUFFC 2009;56(3):489-506
%
% Inputs:
%   rfData          - [samples x channels] raw RF data (CPU or GPU array)
%   rxPos_mm        - [1 x nRx] receive element positions [mm]
%   steerAngle_deg  - Plane wave steering angle [degrees]
%   xGrid           - [1 x nX] lateral pixel positions [mm]
%   zGrid           - [1 x nZ] axial pixel positions [mm]
%   fs_MHz          - Sampling frequency [MHz]
%   c               - Speed of sound [m/s]
%   depthOffset_mm  - VADA depth offset [mm]
%   delayTable      - (optional) precomputed [nPixels x nChannels] sample
%                     index table (gpuArray). Pass [] to compute.
%
% Output:
%   bfImg      - [nZ x nX] complex beamformed image (CPU)
%
% Precomputed delay usage (critical for speed):
%   dt = beamform_planewave_gpu([], rxPos, angle, x, z, fs, c, d, []);
%   for iFrame = 1:N
%       img = beamform_planewave_gpu(rf(:,:,iFrame), rxPos, angle, x, z, fs, c, d, dt);
%   end

persistent gpu_available
if isempty(gpu_available)
    gpu_available = (gpuDeviceCount > 0);
    if gpu_available
        g = gpuDevice;
        fprintf('  [GPU] Using %s (%.1f GB VRAM)\n', g.Name, g.TotalMemory/1e9);
    end
end

[nSamples, nChannels] = size(rfData);
nX = numel(xGrid);
nZ = numel(zGrid);
nPixels = nX * nZ;

c_mm_us = c * 1e-3;  % mm/us
steerAngle_rad = steerAngle_deg * pi / 180;

% Depth offset: the first RF sample corresponds to the round-trip time
% to depthOffset. Since depthOffset is a one-way distance in mm,
% the recording start time is 2 * depthOffset / c.
t0 = 2 * depthOffset_mm / c_mm_us;

% --- Compute delay table if not provided ---
if isempty(delayTable)
    if gpu_available
        xG = gpuArray(single(xGrid(:)));
        zG = gpuArray(single(zGrid(:)));
        rxG = gpuArray(single(rxPos_mm(:)'));
    else
        xG = single(xGrid(:));
        zG = single(zGrid(:));
        rxG = single(rxPos_mm(:)');
    end
    
    [XX, ZZ] = meshgrid(xG, zG);  % [nZ x nX] on GPU
    xx = XX(:);  % [nPixels x 1]
    zz = ZZ(:);
    
    % TX delay: plane wave
    t_tx = (zz * cos(steerAngle_rad) + xx * sin(steerAngle_rad)) / single(c_mm_us);
    
    % RX delay: spherical [nPixels x nChannels]
    t_rx = sqrt((xx - rxG).^2 + zz.^2) / single(c_mm_us);
    
    % Convert to sample indices
    delayTable = (t_tx + t_rx - single(t0)) * single(fs_MHz) + 1;
    
    % If called with empty rfData, just return the delay table
    if isempty(rfData)
        bfImg = delayTable;
        return;
    end
end

% --- IQ demodulation via Hilbert transform ---
if gpu_available && ~isa(rfData, 'gpuArray')
    rfData_gpu = gpuArray(single(rfData));
else
    rfData_gpu = single(rfData);
end

% Hilbert transform on GPU
rfIQ = hilbert_gpu(rfData_gpu);

% --- Vectorized delay-and-sum interpolation on GPU ---
sampleFloor = floor(delayTable);
frac = delayTable - sampleFloor;

% Clamp indices
valid = (sampleFloor >= 1) & (sampleFloor < nSamples);
sampleFloor = max(min(sampleFloor, single(nSamples - 1)), single(1));

% Accumulate across channels
if gpu_available
    pixelVals = gpuArray(zeros(nPixels, 1, 'single'));
else
    pixelVals = zeros(nPixels, 1, 'single');
end

% Vectorized interpolation across all pixels simultaneously per channel
for ich = 1:nChannels
    sf = sampleFloor(:, ich);
    f  = frac(:, ich);
    v  = valid(:, ich);
    
    % Linear interpolation
    val = (1 - f) .* rfIQ(sf, ich) + f .* rfIQ(min(sf + 1, single(nSamples)), ich);
    val(~v) = 0;
    pixelVals = pixelVals + val;
end

% Reshape and gather from GPU
bfImg = reshape(pixelVals, nZ, nX);
if gpu_available
    bfImg = gather(bfImg);
end

end


function xIQ = hilbert_gpu(x)
% HILBERT_GPU  Analytic signal via FFT-based Hilbert transform on GPU.
%
% Equivalent to MATLAB's hilbert() but works on gpuArray.

[n, nCh] = size(x);

% FFT along columns
X = fft(x, [], 1);

% Build multiplier: [2 for positive freq, 0 for negative, 1 for DC and Nyquist]
h = zeros(n, 1, 'single');
if isa(x, 'gpuArray')
    h = gpuArray(h);
end

if mod(n, 2) == 0
    h(1) = 1;           % DC
    h(2:n/2) = 2;       % Positive frequencies
    h(n/2+1) = 1;       % Nyquist
    % h(n/2+2:end) = 0; % Negative frequencies (already zero)
else
    h(1) = 1;
    h(2:(n+1)/2) = 2;
end

xIQ = ifft(X .* h, [], 1);

end
