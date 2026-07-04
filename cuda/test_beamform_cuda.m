%% test_beamform_cuda.m
%  Comprehensive 4-level test suite for the CUDA DAS beamformer.
%
%  Level 1 - Correctness:  Peak locations, cross-validation vs MATLAB
%  Level 2 - Robustness:   Edge cases, bad inputs, boundary conditions
%  Level 3 - Stability:    Memory leaks, reproducibility, long runs
%  Level 4 - Performance:  Benchmarks for both probes
%
%  Run from the cuda/ directory:
%    >> cd cuda
%    >> test_beamform_cuda

clearvars; close all; clc;

%% 0. Compile
fprintf('=== CUDA Beamformer Test Suite ===\n\n');
fprintf('[compile] Building beamform_pw_das.cu ...\n');
try
    mexcuda -R2018a beamform_pw_das.cu -lcufft
    fprintf('[compile] Success.\n\n');
catch ME
    fprintf('[compile] FAILED: %s\n', ME.message);
    fprintf('[compile] If cuFFT not found, try:\n');
    fprintf('  mexcuda -R2018a beamform_pw_das.cu -L"<CUDA_PATH>\\lib\\x64" -lcufft\n\n');
    return;
end

results = struct('name', {}, 'level', {}, 'pass', {}, 'msg', {});

%% ======================================================================
%  LEVEL 1: CORRECTNESS
%  ======================================================================
fprintf('--- Level 1: Correctness ---\n\n');

% 1.1 Point scatterer (L38xp)
[p, m] = test_point_scatterer('L38xp');
results(end+1) = struct('name','1.1 Point scatterer (L38xp)', 'level',1, 'pass',p, 'msg',m);

% 1.2 Point scatterer (UHF29x)
[p, m] = test_point_scatterer('UHF29x');
results(end+1) = struct('name','1.2 Point scatterer (UHF29x)', 'level',1, 'pass',p, 'msg',m);

% 1.3 Multiple scatterers
[p, m] = test_multi_scatterer();
results(end+1) = struct('name','1.3 Multi-scatterer', 'level',1, 'pass',p, 'msg',m);

% 1.4 Steered angle
[p, m] = test_steered_angle();
results(end+1) = struct('name','1.4 Steered angle (+10 deg)', 'level',1, 'pass',p, 'msg',m);

% 1.5 Cross-validation vs MATLAB beamformer
[p, m] = test_cross_validation();
results(end+1) = struct('name','1.5 Cross-validation vs MATLAB', 'level',1, 'pass',p, 'msg',m);

% 1.6 F-number apodization
[p, m] = test_fnum();
results(end+1) = struct('name','1.6 F-number apodization', 'level',1, 'pass',p, 'msg',m);

%% ======================================================================
%  LEVEL 2: ROBUSTNESS
%  ======================================================================
fprintf('\n--- Level 2: Robustness ---\n\n');

% 2.1 Zero input
[p, m] = test_zero_input();
results(end+1) = struct('name','2.1 Zero input', 'level',2, 'pass',p, 'msg',m);

% 2.2 Single frame
[p, m] = test_single_frame();
results(end+1) = struct('name','2.2 Single frame', 'level',2, 'pass',p, 'msg',m);

% 2.3 Multi-frame batch
[p, m] = test_multi_frame();
results(end+1) = struct('name','2.3 Multi-frame (50 frames)', 'level',2, 'pass',p, 'msg',m);

% 2.4 Edge pixel scatterer
[p, m] = test_edge_scatterer();
results(end+1) = struct('name','2.4 Edge pixel scatterer', 'level',2, 'pass',p, 'msg',m);

% 2.5 Bad input dimensions
[p, m] = test_bad_input();
results(end+1) = struct('name','2.5 Bad input (error handling)', 'level',2, 'pass',p, 'msg',m);

%% ======================================================================
%  LEVEL 3: STABILITY (overnight readiness)
%  ======================================================================
fprintf('\n--- Level 3: Stability ---\n\n');

% 3.1 Memory leak test (2000 iterations)
[p, m] = test_memory_leak();
results(end+1) = struct('name','3.1 Memory leak (2000 iters)', 'level',3, 'pass',p, 'msg',m);

% 3.2 Bit-exact reproducibility
[p, m] = test_reproducibility();
results(end+1) = struct('name','3.2 Reproducibility', 'level',3, 'pass',p, 'msg',m);

% 3.3 Long-run correctness (1000 frames)
[p, m] = test_long_run();
results(end+1) = struct('name','3.3 Long-run (1000 frames)', 'level',3, 'pass',p, 'msg',m);

%% ======================================================================
%  LEVEL 4: PERFORMANCE
%  ======================================================================
fprintf('\n--- Level 4: Performance ---\n\n');

% 4.1 Benchmark L38xp
[p, m] = benchmark_probe('L38xp');
results(end+1) = struct('name','4.1 Benchmark L38xp', 'level',4, 'pass',p, 'msg',m);

% 4.2 Benchmark UHF29x
[p, m] = benchmark_probe('UHF29x');
results(end+1) = struct('name','4.2 Benchmark UHF29x', 'level',4, 'pass',p, 'msg',m);

% 4.3 Benchmark vs MATLAB gpuArray
[p, m] = benchmark_vs_matlab();
results(end+1) = struct('name','4.3 CUDA vs MATLAB speedup', 'level',4, 'pass',p, 'msg',m);

%% ======================================================================
%  SUMMARY
%  ======================================================================
fprintf('\n');
fprintf('========================================\n');
fprintf('  TEST SUMMARY\n');
fprintf('========================================\n');
n_pass = sum([results.pass]);
n_total = numel(results);
for i = 1:n_total
    if results(i).pass
        tag = 'PASS';
    else
        tag = 'FAIL';
    end
    fprintf('  [%s] %s -- %s\n', tag, results(i).name, results(i).msg);
end
fprintf('----------------------------------------\n');
fprintf('  %d / %d passed\n', n_pass, n_total);
if n_pass == n_total
    fprintf('  ALL TESTS PASSED\n');
else
    fprintf('  ** %d FAILURE(S) **\n', n_total - n_pass);
end
fprintf('========================================\n');


%% ======================================================================
%  LOCAL FUNCTIONS
%  ======================================================================

% ------ Probe configuration helper ------
function cfg = probe_config(name)
    switch name
        case 'L38xp'
            cfg.nChannels = 64;
            cfg.pitch_mm  = 0.300;
            cfg.fc_hz     = 6e6;
            cfg.fs_hz     = 31.25e6;
            cfg.dx_mm     = 0.050;
            cfg.dz_mm     = 0.025;
            cfg.x_range   = [-8, 8];    % mm
            cfg.z_range   = [2, 18];    % mm
        case 'UHF29x'
            cfg.nChannels = 64;
            cfg.pitch_mm  = 0.090;
            cfg.fc_hz     = 29e6;
            cfg.fs_hz     = 125e6;
            cfg.dx_mm     = 0.011;
            cfg.dz_mm     = 0.005;
            cfg.x_range   = [-2, 2];    % mm
            cfg.z_range   = [1, 5];     % mm
    end
    cfg.name = name;
    cfg.elem_pos_mm = ((0:cfg.nChannels-1) - (cfg.nChannels-1)/2) * cfg.pitch_mm;
end

% ------ Synthetic RF data for point scatterers ------
function [rf, t0_s] = synth_rf(cfg, scatterers, c, tx_angle_rad)
    % scatterers: [N x 2] array of [x_mm, z_mm]
    % Returns rf: [nSamples x nChannels] real single

    if nargin < 4, tx_angle_rad = 0; end

    elem_pos_m = cfg.elem_pos_mm(:) * 1e-3;
    fs = cfg.fs_hz;
    fc = cfg.fc_hz;

    z_min_m = cfg.z_range(1) * 1e-3;
    z_max_m = cfg.z_range(2) * 1e-3;
    t0_s = 2 * z_min_m / c;  % round-trip time to shallowest depth (matches VADA convention)

    % Enough samples to cover deepest round trip
    t_max = z_max_m / c + sqrt(z_max_m^2 + (15e-3)^2) / c;
    nSamples = ceil((t_max - t0_s) * fs) + 200;

    t_axis = t0_s + (0:nSamples-1)' / fs;
    pulse_sigma = 2 / fc;

    rf = zeros(nSamples, cfg.nChannels, 'single');

    for s = 1:size(scatterers, 1)
        xs = scatterers(s, 1) * 1e-3;  % m
        zs = scatterers(s, 2) * 1e-3;  % m

        t_tx = (zs * cos(tx_angle_rad) + xs * sin(tx_angle_rad)) / c;

        for ch = 1:cfg.nChannels
            dx = xs - elem_pos_m(ch);
            dist_rx = sqrt(dx^2 + zs^2);
            t_rx = dist_rx / c;
            t_arrival = t_tx + t_rx;

            dt = t_axis - t_arrival;
            pulse = exp(-dt.^2 / (2 * pulse_sigma^2)) .* cos(2*pi*fc*dt);
            rf(:, ch) = rf(:, ch) + single(pulse);
        end
    end
end

% ------ Convert probe config to MEX args ------
function [elem_m, gx_m, gz_m, fs, t0] = mex_args(cfg, c, t0_s)
    elem_m = single(cfg.elem_pos_mm(:)) * 1e-3;
    gx = cfg.x_range(1):cfg.dx_mm:cfg.x_range(2);
    gz = cfg.z_range(1):cfg.dz_mm:cfg.z_range(2);
    gx_m = single(gx(:)) * 1e-3;
    gz_m = single(gz(:)) * 1e-3;
    fs = cfg.fs_hz;
    t0 = t0_s;
end

% ==================== LEVEL 1 TESTS ====================

function [pass, msg] = test_point_scatterer(probe_name)
    fprintf('  [1] Point scatterer (%s) ... ', probe_name);
    try
        cfg = probe_config(probe_name);
        c = 1540;
        x_scat = 1.0;  z_scat = 8.0;  % mm
        if strcmp(probe_name, 'UHF29x')
            x_scat = 0.5; z_scat = 3.0;
        end

        [rf, t0_s] = synth_rf(cfg, [x_scat, z_scat], c);
        [elem_m, gx_m, gz_m, fs, t0] = mex_args(cfg, c, t0_s);

        bf = beamform_pw_das(rf, elem_m, gx_m, gz_m, 0, c, fs, t0, 0);

        env = abs(bf);
        [~, idx] = max(env(:));
        [iz, ix] = ind2sub(size(env), idx);

        x_peak = gx_m(ix) * 1e3;  % back to mm
        z_peak = gz_m(iz) * 1e3;

        err_x = abs(x_peak - x_scat);
        err_z = abs(z_peak - z_scat);

        tol_x = cfg.dx_mm;
        tol_z = cfg.dz_mm;

        pass = (err_x < tol_x) && (err_z < tol_z);
        msg = sprintf('err=(%.4f, %.4f) mm, tol=(%.4f, %.4f)', ...
            err_x, err_z, tol_x, tol_z);
        if pass, fprintf('PASS\n'); else, fprintf('FAIL: %s\n', msg); end
    catch ME
        pass = false;
        msg = ME.message;
        fprintf('ERROR: %s\n', msg);
    end
end

function [pass, msg] = test_multi_scatterer()
    fprintf('  [1.3] Multi-scatterer ... ');
    try
        cfg = probe_config('L38xp');
        c = 1540;
        scats = [-3, 6; 0, 10; 4, 14; -5, 8; 2, 12];  % 5 scatterers [x_mm, z_mm]

        [rf, t0_s] = synth_rf(cfg, scats, c);
        [elem_m, gx_m, gz_m, fs, t0] = mex_args(cfg, c, t0_s);

        bf = beamform_pw_das(rf, elem_m, gx_m, gz_m, 0, c, fs, t0, 0);
        env = abs(bf);

        % Find peaks: for each scatterer, check local maximum near expected position
        n_found = 0;
        for s = 1:size(scats, 1)
            [~, ix_exp] = min(abs(gx_m*1e3 - scats(s,1)));
            [~, iz_exp] = min(abs(gz_m*1e3 - scats(s,2)));

            % Search in a small ROI around expected position
            roi_x = max(1, ix_exp-10):min(length(gx_m), ix_exp+10);
            roi_z = max(1, iz_exp-10):min(length(gz_m), iz_exp+10);
            roi = env(roi_z, roi_x);

            [~, ri] = max(roi(:));
            [rz, rx] = ind2sub(size(roi), ri);
            iz_peak = roi_z(rz);
            ix_peak = roi_x(rx);

            err_x = abs(gx_m(ix_peak)*1e3 - scats(s,1));
            err_z = abs(gz_m(iz_peak)*1e3 - scats(s,2));

            if err_x < cfg.dx_mm*2 && err_z < cfg.dz_mm*2
                n_found = n_found + 1;
            end
        end

        pass = (n_found == size(scats, 1));
        msg = sprintf('%d/%d scatterers located within 2 pixels', n_found, size(scats,1));
        if pass, fprintf('PASS\n'); else, fprintf('FAIL: %s\n', msg); end
    catch ME
        pass = false; msg = ME.message;
        fprintf('ERROR: %s\n', msg);
    end
end

function [pass, msg] = test_steered_angle()
    fprintf('  [1.4] Steered angle ... ');
    try
        cfg = probe_config('L38xp');
        c = 1540;
        angle_deg = 10;
        angle_rad = angle_deg * pi / 180;
        scats = [0, 10];

        [rf, t0_s] = synth_rf(cfg, scats, c, angle_rad);
        [elem_m, gx_m, gz_m, fs, t0] = mex_args(cfg, c, t0_s);

        bf = beamform_pw_das(rf, elem_m, gx_m, gz_m, angle_rad, c, fs, t0, 0);

        env = abs(bf);
        [~, idx] = max(env(:));
        [iz, ix] = ind2sub(size(env), idx);

        err_x = abs(gx_m(ix)*1e3 - scats(1));
        err_z = abs(gz_m(iz)*1e3 - scats(2));

        pass = (err_x < cfg.dx_mm*2) && (err_z < cfg.dz_mm*2);
        msg = sprintf('angle=%d deg, err=(%.4f, %.4f) mm', angle_deg, err_x, err_z);
        if pass, fprintf('PASS\n'); else, fprintf('FAIL: %s\n', msg); end
    catch ME
        pass = false; msg = ME.message;
        fprintf('ERROR: %s\n', msg);
    end
end

function [pass, msg] = test_cross_validation()
    fprintf('  [1.5] Cross-validation vs MATLAB ... ');
    try
        cfg = probe_config('L38xp');
        c = 1540;
        scats = [2, 10];

        [rf, t0_s] = synth_rf(cfg, scats, c);

        % CUDA path (SI units)
        elem_m = single(cfg.elem_pos_mm(:)) * 1e-3;
        gx_mm = cfg.x_range(1):cfg.dx_mm:cfg.x_range(2);
        gz_mm = cfg.z_range(1):cfg.dz_mm:cfg.z_range(2);
        gx_m = single(gx_mm(:)) * 1e-3;
        gz_m = single(gz_mm(:)) * 1e-3;

        bf_cuda = beamform_pw_das(rf, elem_m, gx_m, gz_m, ...
            0, c, cfg.fs_hz, t0_s, 0);

        % MATLAB path (mm/MHz units)
        depthOffset_mm = cfg.z_range(1);
        bf_matlab = beamform_planewave_gpu(rf, cfg.elem_pos_mm, 0, ...
            gx_mm, gz_mm, cfg.fs_hz/1e6, c, depthOffset_mm, []);

        % Compare peak locations
        env_c = abs(bf_cuda);
        env_m = abs(bf_matlab);

        [~, ic] = max(env_c(:));
        [~, im] = max(env_m(:));
        [izc, ixc] = ind2sub(size(env_c), ic);
        [izm, ixm] = ind2sub(size(env_m), im);

        peak_match = (abs(ixc - ixm) <= 1) && (abs(izc - izm) <= 1);

        % Envelope correlation (normalized)
        ec = env_c(:) / max(env_c(:));
        em = env_m(:) / max(env_m(:));
        corr_val = ec' * em / (norm(ec) * norm(em));

        pass = peak_match && (corr_val > 0.95);
        msg = sprintf('peak offset=(%d,%d) px, envelope corr=%.4f', ...
            abs(ixc-ixm), abs(izc-izm), corr_val);
        if pass, fprintf('PASS\n'); else, fprintf('FAIL: %s\n', msg); end
    catch ME
        pass = false; msg = ME.message;
        fprintf('ERROR: %s\n', msg);
    end
end

function [pass, msg] = test_fnum()
    fprintf('  [1.6] F-number apodization ... ');
    try
        cfg = probe_config('L38xp');
        c = 1540;
        scats = [0, 4];  % shallow scatterer where f-number matters most

        [rf, t0_s] = synth_rf(cfg, scats, c);
        [elem_m, gx_m, gz_m, fs, t0] = mex_args(cfg, c, t0_s);

        % Without f-number
        bf_nofnum = beamform_pw_das(rf, elem_m, gx_m, gz_m, 0, c, fs, t0, 0);
        % With f-number = 1.5
        bf_fnum = beamform_pw_das(rf, elem_m, gx_m, gz_m, 0, c, fs, t0, 1.5);

        env_no = abs(bf_nofnum);
        env_fn = abs(bf_fnum);

        % Both should find the scatterer at the same location
        [~, i1] = max(env_no(:));
        [~, i2] = max(env_fn(:));
        [iz1, ix1] = ind2sub(size(env_no), i1);
        [iz2, ix2] = ind2sub(size(env_fn), i2);

        % F-number reduces active aperture at shallow depths, which can
        % shift the PSF peak by a few pixels. Allow up to 3 pixels.
        peak_match = (abs(ix1-ix2) <= 3) && (abs(iz1-iz2) <= 3);

        % F-number image should have different sidelobe structure
        images_differ = max(abs(env_no(:) - env_fn(:))) > 0.01 * max(env_no(:));

        pass = peak_match && images_differ;
        msg = sprintf('peak offset=(%d,%d) px, images differ=%d', ...
            abs(ix1-ix2), abs(iz1-iz2), images_differ);
        if pass, fprintf('PASS\n'); else, fprintf('FAIL: %s\n', msg); end
    catch ME
        pass = false; msg = ME.message;
        fprintf('ERROR: %s\n', msg);
    end
end

% ==================== LEVEL 2 TESTS ====================

function [pass, msg] = test_zero_input()
    fprintf('  [2.1] Zero input ... ');
    try
        cfg = probe_config('L38xp');
        c = 1540;
        t0_s = cfg.z_range(1) * 1e-3 / c;
        [elem_m, gx_m, gz_m, fs, t0] = mex_args(cfg, c, t0_s);

        nSamples = 500;
        rf = zeros(nSamples, cfg.nChannels, 'single');

        bf = beamform_pw_das(rf, elem_m, gx_m, gz_m, 0, c, fs, t0, 0);

        max_val = max(abs(bf(:)));
        pass = (max_val == 0);
        msg = sprintf('max output = %.2e', max_val);
        if pass, fprintf('PASS\n'); else, fprintf('FAIL: %s\n', msg); end
    catch ME
        pass = false; msg = ME.message;
        fprintf('ERROR: %s\n', msg);
    end
end

function [pass, msg] = test_single_frame()
    fprintf('  [2.2] Single frame ... ');
    try
        cfg = probe_config('L38xp');
        c = 1540;
        [rf, t0_s] = synth_rf(cfg, [0, 10], c);
        [elem_m, gx_m, gz_m, fs, t0] = mex_args(cfg, c, t0_s);

        bf = beamform_pw_das(rf, elem_m, gx_m, gz_m, 0, c, fs, t0, 0);

        pass = (ndims(bf) == 2) && ~isempty(bf) && all(isfinite(bf(:)));
        msg = sprintf('output size = [%s]', num2str(size(bf)));
        if pass, fprintf('PASS\n'); else, fprintf('FAIL: %s\n', msg); end
    catch ME
        pass = false; msg = ME.message;
        fprintf('ERROR: %s\n', msg);
    end
end

function [pass, msg] = test_multi_frame()
    fprintf('  [2.3] Multi-frame ... ');
    try
        cfg = probe_config('L38xp');
        c = 1540;
        [rf1, t0_s] = synth_rf(cfg, [0, 10], c);
        nF = 50;
        rf_batch = repmat(rf1, 1, 1, nF);  % [nSamples x nCh x 50]

        [elem_m, gx_m, gz_m, fs, t0] = mex_args(cfg, c, t0_s);

        bf = beamform_pw_das(rf_batch, elem_m, gx_m, gz_m, 0, c, fs, t0, 0);

        nZ = length(gz_m); nX = length(gx_m);
        correct_size = isequal(size(bf), [nZ, nX, nF]);

        % All frames should be identical (same input)
        frame_diff = max(abs(bf(:,:,1) - bf(:,:,nF)), [], 'all');
        frames_match = (frame_diff == 0);

        pass = correct_size && frames_match;
        msg = sprintf('size=[%dx%dx%d], frame diff=%.2e', ...
            size(bf,1), size(bf,2), size(bf,3), frame_diff);
        if pass, fprintf('PASS\n'); else, fprintf('FAIL: %s\n', msg); end
    catch ME
        pass = false; msg = ME.message;
        fprintf('ERROR: %s\n', msg);
    end
end

function [pass, msg] = test_edge_scatterer()
    fprintf('  [2.4] Edge scatterer ... ');
    try
        cfg = probe_config('L38xp');
        c = 1540;
        % Place scatterer at the edge of the grid
        x_edge = cfg.x_range(2) - cfg.dx_mm;
        z_edge = cfg.z_range(1) + cfg.dz_mm * 5;

        [rf, t0_s] = synth_rf(cfg, [x_edge, z_edge], c);
        [elem_m, gx_m, gz_m, fs, t0] = mex_args(cfg, c, t0_s);

        bf = beamform_pw_das(rf, elem_m, gx_m, gz_m, 0, c, fs, t0, 0);

        pass = all(isfinite(bf(:)));
        msg = sprintf('scatterer at edge (%.1f, %.1f) mm, all finite=%d', ...
            x_edge, z_edge, pass);
        if pass, fprintf('PASS\n'); else, fprintf('FAIL: %s\n', msg); end
    catch ME
        pass = false; msg = ME.message;
        fprintf('ERROR: %s\n', msg);
    end
end

function [pass, msg] = test_bad_input()
    fprintf('  [2.5] Bad input ... ');
    try
        cfg = probe_config('L38xp');
        c = 1540;

        errors_caught = 0;
        expected_errors = 3;

        % Wrong type (double instead of single)
        try
            beamform_pw_das(double(zeros(100,64)), single(zeros(64,1)), ...
                single(zeros(10,1)), single(zeros(10,1)), 0, c, 31.25e6, 0, 0);
        catch
            errors_caught = errors_caught + 1;
        end

        % Mismatched channels
        try
            beamform_pw_das(single(zeros(100,64)), single(zeros(32,1)), ...
                single(zeros(10,1)), single(zeros(10,1)), 0, c, 31.25e6, 0, 0);
        catch
            errors_caught = errors_caught + 1;
        end

        % Too few samples
        try
            beamform_pw_das(single(zeros(2,64)), single(zeros(64,1)), ...
                single(zeros(10,1)), single(zeros(10,1)), 0, c, 31.25e6, 0, 0);
        catch
            errors_caught = errors_caught + 1;
        end

        pass = (errors_caught == expected_errors);
        msg = sprintf('%d/%d bad inputs caught cleanly', errors_caught, expected_errors);
        if pass, fprintf('PASS\n'); else, fprintf('FAIL: %s\n', msg); end
    catch ME
        pass = false; msg = ME.message;
        fprintf('ERROR: %s\n', msg);
    end
end

% ==================== LEVEL 3 TESTS ====================

function [pass, msg] = test_memory_leak()
    fprintf('  [3.1] Memory leak test (2000 iters) ... ');
    try
        cfg = probe_config('L38xp');
        c = 1540;

        % Small data for fast iterations
        nSamples = 500;
        nCh = 64;
        rf = randn(nSamples, nCh, 'single');
        elem_m = single(cfg.elem_pos_mm(:)) * 1e-3;
        gx_m = single(linspace(-5e-3, 5e-3, 80)');
        gz_m = single(linspace(2e-3, 10e-3, 80)');
        t0 = 2e-3 / c;

        % Warmup (first call may allocate cuFFT caches)
        for i = 1:5
            bf = beamform_pw_das(rf, elem_m, gx_m, gz_m, 0, c, cfg.fs_hz, t0, 0);
        end
        clear bf;
        g = gpuDevice;
        wait(g);
        mem_before = g.FreeMemory;

        % Run 2000 iterations
        nIter = 2000;
        for i = 1:nIter
            bf = beamform_pw_das(rf, elem_m, gx_m, gz_m, 0, c, cfg.fs_hz, t0, 0);
        end
        clear bf;
        wait(g);
        mem_after = g.FreeMemory;

        leak_mb = (double(mem_before) - double(mem_after)) / 1e6;
        pass = abs(leak_mb) < 2.0;  % tolerance: 2 MB
        msg = sprintf('drift = %.2f MB over %d iters', leak_mb, nIter);
        if pass, fprintf('PASS\n'); else, fprintf('FAIL: %s\n', msg); end
    catch ME
        pass = false; msg = ME.message;
        fprintf('ERROR: %s\n', msg);
    end
end

function [pass, msg] = test_reproducibility()
    fprintf('  [3.2] Reproducibility ... ');
    try
        cfg = probe_config('L38xp');
        c = 1540;
        [rf, t0_s] = synth_rf(cfg, [1, 10], c);
        [elem_m, gx_m, gz_m, fs, t0] = mex_args(cfg, c, t0_s);

        bf1 = beamform_pw_das(rf, elem_m, gx_m, gz_m, 0, c, fs, t0, 0);
        bf2 = beamform_pw_das(rf, elem_m, gx_m, gz_m, 0, c, fs, t0, 0);

        max_diff = max(abs(bf1(:) - bf2(:)));
        pass = (max_diff == 0);
        msg = sprintf('max diff between runs = %.2e', max_diff);
        if pass, fprintf('PASS\n'); else, fprintf('FAIL: %s\n', msg); end
    catch ME
        pass = false; msg = ME.message;
        fprintf('ERROR: %s\n', msg);
    end
end

function [pass, msg] = test_long_run()
    fprintf('  [3.3] Long-run correctness (1000 frames) ... ');
    try
        cfg = probe_config('L38xp');
        c = 1540;

        nSamples = 500;
        rf_single = randn(nSamples, cfg.nChannels, 'single');
        rf_batch = repmat(rf_single, 1, 1, 1000);

        elem_m = single(cfg.elem_pos_mm(:)) * 1e-3;
        gx_m = single(linspace(-5e-3, 5e-3, 60)');
        gz_m = single(linspace(2e-3, 10e-3, 60)');
        t0 = 2e-3 / c;

        bf = beamform_pw_das(rf_batch, elem_m, gx_m, gz_m, ...
            0, c, cfg.fs_hz, t0, 0);

        % Compare first and last frame
        diff_first_last = max(abs(bf(:,:,1) - bf(:,:,1000)), [], 'all');

        % Compare first and middle frame
        diff_first_mid = max(abs(bf(:,:,1) - bf(:,:,500)), [], 'all');

        pass = (diff_first_last == 0) && (diff_first_mid == 0);
        msg = sprintf('frame 1 vs 500 diff=%.2e, vs 1000 diff=%.2e', ...
            diff_first_mid, diff_first_last);
        if pass, fprintf('PASS\n'); else, fprintf('FAIL: %s\n', msg); end
    catch ME
        pass = false; msg = ME.message;
        fprintf('ERROR: %s\n', msg);
    end
end

% ==================== LEVEL 4 BENCHMARKS ====================

function [pass, msg] = benchmark_probe(probe_name)
    fprintf('  [4] Benchmark %s ... ', probe_name);
    try
        cfg = probe_config(probe_name);
        c = 1540;

        [rf1, t0_s] = synth_rf(cfg, [0, 8], c);
        nSamples = size(rf1, 1);
        nF = 100;
        rf_batch = repmat(rf1, 1, 1, nF);

        [elem_m, gx_m, gz_m, fs, t0] = mex_args(cfg, c, t0_s);
        nZ = length(gz_m); nX = length(gx_m);

        % Warmup
        beamform_pw_das(rf1, elem_m, gx_m, gz_m, 0, c, fs, t0, 0);

        % Benchmark
        tic;
        bf = beamform_pw_das(rf_batch, elem_m, gx_m, gz_m, 0, c, fs, t0, 0);
        t = toc;

        fps = nF / t;
        mpx_s = nZ * nX * nF / t / 1e6;

        pass = true;  % benchmarks always "pass"
        msg = sprintf('%d frames in %.2fs (%.0f fps, %.1f Mpx/s, grid %dx%d)', ...
            nF, t, fps, mpx_s, nZ, nX);
        fprintf('%s\n', msg);
    catch ME
        pass = false; msg = ME.message;
        fprintf('ERROR: %s\n', msg);
    end
end

function [pass, msg] = benchmark_vs_matlab()
    fprintf('  [4.3] CUDA vs MATLAB speedup ... ');
    try
        cfg = probe_config('L38xp');
        c = 1540;
        [rf, t0_s] = synth_rf(cfg, [0, 10], c);

        % CUDA timing
        elem_m = single(cfg.elem_pos_mm(:)) * 1e-3;
        gx_mm = cfg.x_range(1):cfg.dx_mm:cfg.x_range(2);
        gz_mm = cfg.z_range(1):cfg.dz_mm:cfg.z_range(2);
        gx_m = single(gx_mm(:)) * 1e-3;
        gz_m = single(gz_mm(:)) * 1e-3;
        depthOffset_mm = cfg.z_range(1);

        % Warmup both
        beamform_pw_das(rf, elem_m, gx_m, gz_m, 0, c, cfg.fs_hz, t0_s, 0);
        beamform_planewave_gpu(rf, cfg.elem_pos_mm, 0, gx_mm, gz_mm, ...
            cfg.fs_hz/1e6, c, depthOffset_mm, []);

        nRuns = 20;

        tic;
        for i = 1:nRuns
            bf_c = beamform_pw_das(rf, elem_m, gx_m, gz_m, 0, c, cfg.fs_hz, t0_s, 0);
        end
        t_cuda = toc / nRuns;

        tic;
        for i = 1:nRuns
            bf_m = beamform_planewave_gpu(rf, cfg.elem_pos_mm, 0, gx_mm, gz_mm, ...
                cfg.fs_hz/1e6, c, depthOffset_mm, []);
        end
        t_matlab = toc / nRuns;

        speedup = t_matlab / t_cuda;

        pass = true;
        msg = sprintf('CUDA=%.1fms, MATLAB=%.1fms, speedup=%.1fx', ...
            t_cuda*1e3, t_matlab*1e3, speedup);
        fprintf('%s\n', msg);
    catch ME
        pass = false; msg = ME.message;
        fprintf('ERROR: %s\n', msg);
    end
end
