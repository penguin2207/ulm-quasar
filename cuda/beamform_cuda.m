function bfImg = beamform_cuda(rfData, rxPos_mm, steerAngle_deg, ...
    xGrid, zGrid, fs_MHz, c, depthOffset_mm, delayTable, fnum)
% BEAMFORM_CUDA  CUDA-accelerated delay-and-sum beamforming.
%
% Drop-in replacement for beamform_planewave_gpu.  Uses a MEX-CUDA kernel
% with Catmull-Rom cubic interpolation, built-in Hilbert transform, and
% optional f-number apodization.
%
% Signature matches beamform_planewave_gpu (9 args) plus optional fnum:
%
%   bfImg = beamform_cuda(rfData, rxPos_mm, steerAngle_deg, ...
%               xGrid, zGrid, fs_MHz, c, depthOffset_mm, delayTable)
%
%   bfImg = beamform_cuda(..., fnum)      % with f-number apodization
%
% The delayTable argument is accepted for API compatibility but ignored
% by the CUDA path (delays are computed on-the-fly per thread).
%
% When rfData is empty, computes and returns a delay table via the MATLAB
% gpuArray path (called once at pipeline startup; not performance-critical).
%
% Inputs:
%   rfData          [nSamples x nChannels] real single (raw RF)
%   rxPos_mm        [1 x nRx] receive element positions [mm]
%   steerAngle_deg  scalar, steering angle [degrees]
%   xGrid           [1 x nX] lateral pixel positions [mm]
%   zGrid           [1 x nZ] axial pixel positions [mm]
%   fs_MHz          scalar, sampling frequency [MHz]
%   c               scalar, speed of sound [m/s]
%   depthOffset_mm  scalar, VADA depth offset [mm]
%   delayTable      precomputed delay table (ignored by CUDA path)
%   fnum            (optional) f-number for apodization, <=0 to disable
%
% Output:
%   bfImg           [nZ x nX] complex single (beamformed IQ, on CPU)

% Default: no f-number apodization (matches legacy behavior)
if nargin < 10 || isempty(fnum)
    fnum = 0;
end

% ---- Check MEX availability (once per session) ----
persistent mex_ok
if isempty(mex_ok)
    mex_ok = (exist('beamform_pw_das', 'file') == 3);
    if mex_ok
        fprintf('  [CUDA] MEX beamformer loaded\n');
    else
        fprintf('  [CUDA] MEX not found -- falling back to gpuArray beamformer\n');
    end
end

% ---- Delay table precomputation (empty rfData) ----
% This path is called once at pipeline startup. The returned table is
% passed back on subsequent calls but ignored by the CUDA kernel (it
% computes delays on-the-fly). We delegate to the MATLAB beamformer
% so the existing pipeline code works without changes.
if isempty(rfData)
    bfImg = beamform_planewave_gpu([], rxPos_mm, steerAngle_deg, ...
        xGrid, zGrid, fs_MHz, c, depthOffset_mm, []);
    return;
end

% ---- Fallback to MATLAB if MEX not compiled ----
if ~mex_ok
    bfImg = beamform_planewave_gpu(rfData, rxPos_mm, steerAngle_deg, ...
        xGrid, zGrid, fs_MHz, c, depthOffset_mm, delayTable);
    return;
end

% ---- Unit conversion to SI (MEX expects meters, Hz, seconds) ----
elem_m    = single(rxPos_mm(:)) * 1e-3;            % mm -> m
gx_m      = single(xGrid(:))   * 1e-3;             % mm -> m
gz_m      = single(zGrid(:))   * 1e-3;             % mm -> m
angle_rad = double(steerAngle_deg) * pi / 180;     % deg -> rad
c_mps     = double(c);                              % already m/s
fs_hz     = double(fs_MHz) * 1e6;                   % MHz -> Hz
t0_s      = 2 * double(depthOffset_mm) * 1e-3 / c_mps;  % round-trip start time [s]

% Ensure input is real single on CPU (gather handles gpuArray)
rf = single(gather(real(rfData)));

% ---- Call CUDA MEX ----
bfImg = beamform_pw_das(rf, elem_m, gx_m, gz_m, ...
    angle_rad, c_mps, fs_hz, t0_s, double(fnum));

end
