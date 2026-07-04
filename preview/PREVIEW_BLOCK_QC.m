function qcTable = PREVIEW_BLOCK_QC(folder, pitch_mm, xRange, zRange, filterPattern)
%PREVIEW_BLOCK_QC  Fast per-block signal-level QC for a folder of VADA.
%
%  Usage:
%    qc = PREVIEW_BLOCK_QC(folder, pitch_mm, xRange, zRange)              -- all .vada
%    qc = PREVIEW_BLOCK_QC(folder, pitch_mm, xRange, zRange, 'C1b')       -- filter to 'C1b'
%
%  For each .vada in `folder` (optionally name-filtered by `filterPattern`),
%  beamform a small (200-frame) ensemble and compute: peak power, mean
%  power, pre-SVD-3 support density. Returns a table sorted by peak power.
%  Useful for identifying a "dim" block that might be an unlabeled background.
%
%  Output qcTable columns:
%    basename, peakPower_dB, meanPower_dB, supportDensity_pct, probableBg
%
%  `probableBg = true` when peakPower is >= 6 dB below median of siblings.

if nargin < 4, error('Usage: PREVIEW_BLOCK_QC(folder, pitch_mm, xRange, zRange [, filterPattern])'); end
if nargin < 5, filterPattern = ''; end

addpath('cuda');
vadaScriptsPath = 'C:\path\to\VevoF2\RFExampleScriptsF2\Matlab';
if exist(vadaScriptsPath, 'dir'), addpath(genpath(vadaScriptsPath)); end

modeName = '.vada';
c        = 1540;
N_SAMPLE = 200;              % 200-frame sample per block
SVD_CUT  = 3;
USE_GPU  = true;

vf = dir(fullfile(folder, '*.vada'));
vf = vf(~contains({vf.name}, 'bg', 'IgnoreCase', true));
if ~isempty(filterPattern)
    vf = vf(contains({vf.name}, filterPattern) & ~contains({vf.name}, '2x'));
end
if isempty(vf), error('No matching .vada files in %s (filter=''%s'')', folder, filterPattern); end

dx = 0.040; dz = 0.040;
if pitch_mm < 0.2, dx = 0.020; dz = 0.020; end  % UHF29x finer grid
xGrid = xRange(1):dx:xRange(2);
zGrid = zRange(1):dz:zRange(2);
nX = numel(xGrid); nZ = numel(zGrid);

N = numel(vf);
peakPwr = zeros(N,1);
meanPwr = zeros(N,1);
supDensity = zeros(N,1);
names = cell(N,1);

fprintf('\n=== PREVIEW_BLOCK_QC (%d blocks) ===\n', N);
for i = 1:N
    basename = regexprep(vf(i).name, '\.vada$', '');
    names{i} = basename;
    try
        meta = acq_load_block_meta(folder, basename, modeName, c, pitch_mm);
        nFrames = min(N_SAMPLE, meta.numCompoundFrames);

        rxPos_mm  = meta.elemPos_mm(meta.anglePairs(meta.zeroAngleIdx).rxElements);
        angle_deg = meta.anglePairs(meta.zeroAngleIdx).angle;

        evList = (0:nFrames-1)*meta.eventsPerFrame + meta.zeroEventIdx;
        [VadaChunk, ~, ~, ~] = VsiVadaDataRead(folder, basename, evList, modeName);
        nSrf = size(VadaChunk(1).Data,1); nRrf = size(VadaChunk(1).Data,2);
        rfBatch = zeros(nSrf, nRrf, nFrames, 'single');
        for k = 1:nFrames, rfBatch(:,:,k) = single(VadaChunk(k).Data); end

        IQ = beamform_cuda(rfBatch, rxPos_mm, angle_deg, xGrid, zGrid, ...
            meta.fs_MHz, meta.c, meta.depthOffset_mm, []);
        IQ = single(IQ);
        clear rfBatch VadaChunk;

        IQf = svd_clutter_filter_rsvd(IQ, SVD_CUT, [], USE_GPU);
        pwr = mean(abs(IQf).^2, 3);
        clear IQ IQf;

        p_norm = pwr / max(pwr(:) + eps);
        p_db = 10*log10(p_norm + 1e-6);
        peakPwr(i) = 10*log10(max(pwr(:)) + eps);
        meanPwr(i) = 10*log10(mean(pwr(:)) + eps);

        thr = 0.1 * max(pwr(:));
        supDensity(i) = 100 * nnz(pwr > thr) / numel(pwr);

        fprintf('  [%2d/%2d] %s  peak=%.1f dB  mean=%.1f dB  sup@10%%=%.1f%%\n', ...
            i, N, basename, peakPwr(i), meanPwr(i), supDensity(i));
    catch ME
        fprintf('  [%2d/%2d] %s  ERROR: %s\n', i, N, basename, ME.message);
        peakPwr(i) = NaN; meanPwr(i) = NaN; supDensity(i) = NaN;
    end
end

probableBg = (peakPwr < nanmedian(peakPwr) - 6);

qcTable = table(names, peakPwr, meanPwr, supDensity, probableBg, ...
    'VariableNames', {'basename', 'peakPower_dB', 'meanPower_dB', 'supportDensity_pct', 'probableBg'});
qcTable = sortrows(qcTable, 'peakPower_dB');

fprintf('\n--- Ranked by peak power (lowest first) ---\n');
disp(qcTable);

probableList = qcTable.basename(qcTable.probableBg);
if ~isempty(probableList)
    fprintf('\n*** PROBABLE BACKGROUND (>=6 dB below sibling median): ***\n');
    for i = 1:numel(probableList)
        fprintf('  %s\n', probableList{i});
    end
else
    fprintf('\nAll blocks within 6 dB of sibling median -- no probable-bg outlier.\n');
end
end
