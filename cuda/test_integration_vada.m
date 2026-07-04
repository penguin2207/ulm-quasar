%% test_integration_vada.m
%  Integration test: compare CUDA and MATLAB beamformers on real VADA data.
%
%  This is the final gate before overnight use. It loads one block of
%  real acquisition data and runs both beamformers side by side.
%
%  BEFORE RUNNING: set dataFolder and baseFilename below to point at
%  a real VADA dataset.

clearvars; close all; clc;
fprintf('=== Integration Test: CUDA vs MATLAB on Real VADA Data ===\n\n');

%% ---- Configuration ----
% CHANGE THESE to point at a real dataset:
dataFolder   = '';   % e.g. 'D:\Data\2026-04-08\...'
baseFilename = '';   % e.g. 'MB_Ladder_001'
modeName     = 'raw';

if isempty(dataFolder) || isempty(baseFilename)
    fprintf('ERROR: Set dataFolder and baseFilename before running.\n');
    fprintf('       Edit the top of test_integration_vada.m.\n');
    return;
end

%% ---- Load first compound frame ----
fprintf('[load] Reading VADA metadata ...\n');
numEvents = 6;  % typical: 3 angles x 2 polarities, or 1 angle x 2 pol
[VadaTest, Param, TxrParam, ~] = VsiVadaDataRead(dataFolder, baseFilename, ...
    1:numEvents, modeName);

% Auto-detect probe parameters
if TxrParam.ArrayPitch == 0
    pitch_mm = 0.300;  % fallback for L38xp
    fprintf('  WARNING: ArrayPitch=0 in metadata, using %.3f mm\n', pitch_mm);
elseif TxrParam.ArrayPitch >= 10
    pitch_mm = TxrParam.ArrayPitch / 1000;
else
    pitch_mm = TxrParam.ArrayPitch;
end

c = 1540;  % m/s (override for agarose/tissue)
fs_MHz = Param.SampleFreq;
depthOffset_mm = Param.DepthOffset_mm;
nRx = size(VadaTest(1).Data, 2);

fprintf('  Probe: %s | %d RX elem | pitch=%.4f mm | Fs=%.1f MHz\n', ...
    TxrParam.Name, nRx, pitch_mm, fs_MHz);
fprintf('  Depth offset: %.2f mm | SoS: %d m/s\n', depthOffset_mm, c);

% Element positions (centered)
elemPos_mm = ((1:nRx) - (nRx+1)/2) * pitch_mm;

% Beamforming grid (auto-scaled from wavelength)
txFreq_MHz = VadaTest(1).Waveform.Channel(1).frequency;
lambda_mm = c * 1e-3 / txFreq_MHz;
dx_mm = round(lambda_mm / 5, 4);
dz_mm = round(lambda_mm / 10, 4);

rxSpan = (nRx - 1) * pitch_mm;
xRange = [min(elemPos_mm) - rxSpan*0.2, max(elemPos_mm) + rxSpan*0.2];
zRange_mm = [depthOffset_mm, depthOffset_mm + ...
    size(VadaTest(1).Data, 1) / (fs_MHz * 2) * c * 1e-3];

xGrid = xRange(1):dx_mm:xRange(2);
zGrid = zRange_mm(1):dz_mm:zRange_mm(2);

fprintf('  Grid: %d x %d pixels (dx=%.4f mm, dz=%.4f mm)\n', ...
    length(zGrid), length(xGrid), dx_mm, dz_mm);

% Use zero-angle event (or first event if single-angle)
rfData = single(VadaTest(1).Data);  % [nSamples x nRx]
fprintf('  RF data: [%d x %d] samples\n\n', size(rfData));

%% ---- Run MATLAB beamformer ----
fprintf('[matlab] Running beamform_planewave_gpu ...\n');
tic;
bf_matlab = beamform_planewave_gpu(rfData, elemPos_mm, 0, ...
    xGrid, zGrid, fs_MHz, c, depthOffset_mm, []);
t_matlab = toc;
fprintf('  Time: %.3f s\n', t_matlab);

%% ---- Run CUDA beamformer ----
fprintf('[cuda] Running beamform_pw_das ...\n');

% Convert to SI units
elem_m = single(elemPos_mm(:)) * 1e-3;
gx_m   = single(xGrid(:)) * 1e-3;
gz_m   = single(zGrid(:)) * 1e-3;
t0_s   = 2 * depthOffset_mm * 1e-3 / c;

tic;
bf_cuda = beamform_pw_das(rfData, elem_m, gx_m, gz_m, ...
    0, double(c), double(fs_MHz * 1e6), t0_s, 0);
t_cuda = toc;
fprintf('  Time: %.3f s\n', t_cuda);
fprintf('  Speedup: %.1fx\n\n', t_matlab / t_cuda);

%% ---- Compare results ----
fprintf('[compare] Analyzing differences ...\n');

env_c = abs(bf_cuda);
env_m = abs(bf_matlab);

% Peak locations
[~, ic] = max(env_c(:));
[~, im] = max(env_m(:));
[izc, ixc] = ind2sub(size(env_c), ic);
[izm, ixm] = ind2sub(size(env_m), im);
fprintf('  CUDA peak:   (%d, %d) = (%.2f mm, %.2f mm)\n', ...
    izc, ixc, xGrid(ixc), zGrid(izc));
fprintf('  MATLAB peak: (%d, %d) = (%.2f mm, %.2f mm)\n', ...
    izm, ixm, xGrid(ixm), zGrid(izm));
fprintf('  Peak offset: (%d, %d) pixels\n', abs(izc-izm), abs(ixc-ixm));

% Envelope correlation
ec = env_c(:) / max(env_c(:));
em = env_m(:) / max(env_m(:));
corr_val = ec' * em / (norm(ec) * norm(em));
fprintf('  Envelope correlation: %.6f\n', corr_val);

% Max absolute error (normalized)
max_err = max(abs(env_c(:) - env_m(:))) / max(env_m(:));
fprintf('  Max norm. error: %.4f\n', max_err);

% SSIM (structural similarity)
env_c_norm = env_c / max(env_c(:));
env_m_norm = env_m / max(env_m(:));
ssim_val = ssim(env_c_norm, env_m_norm);
fprintf('  SSIM: %.6f\n\n', ssim_val);

%% ---- Visual comparison ----
fprintf('[plot] Generating comparison figures ...\n');

db_range = [-40 0];

figure('Name', 'Integration Test: CUDA vs MATLAB', ...
    'Position', [100 100 1400 500]);

subplot(1,3,1);
env_db_m = 20*log10(env_m / max(env_m(:)));
imagesc(xGrid, zGrid, env_db_m, db_range);
colormap hot; colorbar;
title(sprintf('MATLAB (%.1f ms)', t_matlab*1e3));
xlabel('Lateral (mm)'); ylabel('Depth (mm)');
axis image;

subplot(1,3,2);
env_db_c = 20*log10(env_c / max(env_c(:)));
imagesc(xGrid, zGrid, env_db_c, db_range);
colormap hot; colorbar;
title(sprintf('CUDA (%.1f ms)', t_cuda*1e3));
xlabel('Lateral (mm)'); ylabel('Depth (mm)');
axis image;

subplot(1,3,3);
diff_db = 20*log10(abs(env_c - env_m) / max(env_m(:)) + 1e-10);
imagesc(xGrid, zGrid, diff_db, [-60 -20]);
colormap hot; colorbar;
title(sprintf('Difference (corr=%.4f)', corr_val));
xlabel('Lateral (mm)'); ylabel('Depth (mm)');
axis image;

sgtitle('CUDA vs MATLAB Beamformer Comparison');

%% ---- Pass/Fail ----
fprintf('\n');
peak_ok = (abs(izc-izm) <= 2) && (abs(ixc-ixm) <= 2);
corr_ok = corr_val > 0.95;
ssim_ok = ssim_val > 0.90;

if peak_ok && corr_ok && ssim_ok
    fprintf('RESULT: PASS -- CUDA beamformer matches MATLAB reference.\n');
    fprintf('        Safe for overnight use.\n');
else
    fprintf('RESULT: FAIL --\n');
    if ~peak_ok, fprintf('  Peak location mismatch > 2 pixels\n'); end
    if ~corr_ok, fprintf('  Envelope correlation < 0.95\n'); end
    if ~ssim_ok, fprintf('  SSIM < 0.90\n'); end
    fprintf('        Investigate before overnight use.\n');
end
