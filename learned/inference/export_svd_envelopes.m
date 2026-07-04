% SPDX-License-Identifier: MIT
function export_svd_envelopes(polCachePath, outMat, domain, nFr, useGPU)
% EXPORT_SVD_ENVELOPES  SVD-filtered envelope frames for Python learned-localizer inference.
%
% NEW SPEC (REANALYSIS_INPUTS_RESPONSE.md Sections 2-4): the input is a per-block POLARITY
% cache '<label>_POL.mat' (top-level vars IQ_pos, IQ_neg, nFrames, g, blockMeta). A `domain`
% is derived on the fly (PI / fundamental / singlepol), then run through the IDENTICAL
% reanalysis front end (rSVD, cutoff 8, seed 12345, GPU) so every localizer (classical +
% deconv + Deep-ULM + LISTA) sees the same filtered input the classical anchor saw.
%
%   polCachePath : full path to '<label>_POL.mat' with top-level vars:
%                    IQ_pos  complex single [g.nZ x g.nX x nFrames]
%                    IQ_neg  complex single [...]  OR  []  (empty when ~blockMeta.hasPI)
%                    nFrames scalar double
%                    g       per-block natural-depth grid struct (Section 4a)
%                    blockMeta  compact meta struct
%   outMat       : output .mat (-v7.3) with vars:
%                    env     [g.nZ x g.nX x nFr] single = single(abs(IQf))  (net input)
%                    g       the per-block grid (so downstream uses S.g, NOT a global config)
%                    nFrames the exported frame count (= size(env,3))
%                    domain  the domain string (bookkeeping)
%   domain       : 'PI' (default) | 'fundamental' | 'singlepol'   (Section 2 step 1)
%   nFr          : number of frames to export ([]/missing -> all)
%   useGPU       : rSVD on GPU (default false). LAB-PC FINAL RUN: TRUE to byte-match the
%                  reanalysis envelope (GPU vs CPU changes the rSVD random projection slightly
%                  even with the same seed).
%
% Uses the reanalysis core: core/svd_clutter_filter_rsvd.m (the ulm-pipeline repo core/).

CORE = fullfile(fileparts(mfilename('fullpath')), '..', '..', 'core');  % learned/inference -> repo core/
addpath(CORE);
if nargin < 3 || isempty(domain), domain = 'PI'; end
if nargin < 5 || isempty(useGPU), useGPU = false; end

% ---- load the polarity cache (top-level vars; saved -struct, RESPONSE.md Section 3) ----
S = load(polCachePath);
assert(isfield(S, 'IQ_pos'), 'export_svd_envelopes:noIQpos', ...
    '%s has no top-level IQ_pos', polCachePath);
IQ_pos = S.IQ_pos;
if isfield(S, 'IQ_neg'), IQ_neg = S.IQ_neg; else, IQ_neg = []; end   % handle missing / empty

% ---- derive the domain image stack D (Section 2 step 1; handle isempty(IQ_neg)) ----
switch lower(domain)
    case 'pi'
        if isempty(IQ_neg), D = IQ_pos; else, D = IQ_pos + IQ_neg; end
    case 'fundamental'
        if isempty(IQ_neg), D = IQ_pos; else, D = IQ_pos - IQ_neg; end
    case 'singlepol'
        D = IQ_pos;
    otherwise
        error('export_svd_envelopes:badDomain', ...
            'unknown domain %s (PI | fundamental | singlepol)', domain);
end

% ---- frame cap ----
nTot = size(D, 3);
if nargin < 4 || isempty(nFr), nFr = nTot; end
nFr = min(nFr, nTot);
D = D(:, :, 1:nFr);

% ---- IDENTICAL reanalysis SVD clutter filter: cutoff 8, upper [], seed 12345 (Section 2 step 2)
IQf = svd_clutter_filter_rsvd(D, 8, [], useGPU, 12345);
env = single(abs(IQf));                  % [g.nZ x g.nX x nFr]  -> net / clip-high input

% ---- carry the per-block grid + frame count so downstream uses S.g (Section 4a) ----
g       = S.g;
nFrames = nFr;

save(outMat, 'env', 'g', 'nFrames', 'domain', '-v7.3');
fprintf('[export] %s : env %d x %d x %d  domain=%s  useGPU=%d\n', ...
        outMat, size(env,1), size(env,2), size(env,3), domain, useGPU);
end
