function IQ_filtered = svd_clutter_filter_rsvd(IQ_stack, cutoffLow, cutoffHigh, useGPU, seed)
% SVD_CLUTTER_FILTER_RSVD  Fast randomized-SVD clutter filter.
%
% Drop-in replacement for svd_clutter_filter_gpu using randomized SVD
% (Halko, Martinsson, Tropp 2011) to compute only the top-k singular
% vectors needed for clutter rejection. Typically 5-10x faster than
% full economy SVD when cutoffLow << nFrames.
%
% The speedup comes from a simple observation: the original filter
% computes ALL nFrames singular vectors via svd(S,'econ') but only
% needs the first cutoffLow of them. Randomized SVD computes only
% the top-(cutoffLow + oversample) vectors at a fraction of the cost.
%
% Clutter subtraction uses the projection identity:
%   S_filtered = S - S * V_clutter * V_clutter'
% where V_clutter is the top-cutoffLow right singular vectors.
% This avoids materializing U and Sigma entirely.
%
% Reference: Halko, Martinsson, Tropp. "Finding structure with
% randomness" SIAM Review 53:2 (2011).
%
% Inputs:
%   IQ_stack   - [nZ x nX x nFrames] complex IQ image stack
%   cutoffLow  - Remove top-cutoffLow singular vectors (tissue)
%   cutoffHigh - Remove above cutoffHigh (noise). [] = keep all below.
%                Not supported by the fast path (falls back to full SVD
%                if set, since randomized SVD only yields the top-k).
%   useGPU     - Boolean: use GPU acceleration
%   seed       - (optional) RNG seed for the randomized range finder.
%                [] or omitted = unseeded (nondeterministic, as published).
%                When set, seeds BOTH the CPU (rng) and GPU (gpurng) streams,
%                so the result is bit-reproducible on the same host/session.
%
% Output:
%   IQ_filtered - [nZ x nX x nFrames] filtered stack

[nZ, nX, nFrames] = size(IQ_stack);
nPixels = nZ * nX;

tSVD = tic;

% --- Fall back to full SVD if cutoffHigh is set ---
% Randomized SVD only gives the TOP singular vectors, not the tail.
% The noise-removal case (cutoffHigh set) needs the tail, so the
% fast path isn't applicable. In practice cutoffHigh is almost
% always empty for ULM clutter removal.
if ~isempty(cutoffHigh) && cutoffHigh > 0
    IQ_filtered = svd_clutter_filter_gpu(IQ_stack, cutoffLow, cutoffHigh, useGPU);
    return;
end

% --- Reshape to Casorati matrix [pixels x frames] ---
S = reshape(IQ_stack, nPixels, nFrames);

if useGPU
    S = gpuArray(single(S));
else
    S = single(S);
end

% --- Randomized range finder ---
% Oversampling p=10 is standard (Halko et al.) and gives excellent
% accuracy for the top-cutoffLow subspace. Cost dominates at
% k=cutoffLow+oversample which is still << nFrames for typical ULM.
oversample = 10;
k = min(cutoffLow + oversample, nFrames);

% --- Deterministic seeding (optional) ---
% Seed BOTH streams: rng() controls only the CPU Mersenne-Twister, but when
% useGPU the Omega draw uses 'gpuArray' randn, which is governed by the GPU
% stream and stays nondeterministic unless gpurng() is also set. Seeding only
% rng() would silently leave the GPU path nondeterministic (Tier-A repro fails).
if nargin >= 5 && ~isempty(seed)
    rng(seed);
    if useGPU
        gpurng(seed);
    end
end

% Random test matrix (complex, matching S)
if useGPU
    Omega = complex( ...
        randn(nFrames, k, 'single', 'gpuArray'), ...
        randn(nFrames, k, 'single', 'gpuArray'));
else
    Omega = complex( ...
        randn(nFrames, k, 'single'), ...
        randn(nFrames, k, 'single'));
end

% Sample the range: Y spans the top-k row space of S
Y = S * Omega;          % [nPixels x k]
[Q, ~] = qr(Y, 0);      % [nPixels x k], orthonormal columns

% Project S onto the range -> small k x nFrames matrix
B = Q' * S;             % [k x nFrames]

% SVD of small matrix (trivial cost: k x nFrames is ~15 x 391)
[~, ~, V] = svd(B, 'econ');  % V: [nFrames x k]

% --- Clutter subtraction via projection identity ---
%
% Mathematical identity:
%   S = U*Sigma*V'         (full SVD)
%   Top-k component:      U(:,1:k)*Sigma(1:k,1:k)*V(:,1:k)'
%                       = S * V(:,1:k) * V(:,1:k)'
%   (because  S * V(:,1:k) = U(:,1:k)*Sigma(1:k,1:k))
%
% So to subtract the top-cutoffLow clutter components:
%   S_filtered = S - S * V(:,1:cutoffLow) * V(:,1:cutoffLow)'
%
% Only two small matrix multiplies needed, no U or Sigma required.
V_clutter = V(:, 1:cutoffLow);              % [nFrames x cutoffLow]
SVc = S * V_clutter;                         % [nPixels x cutoffLow]
S_filtered = S - SVc * V_clutter';           % [nPixels x nFrames]

% --- Reshape and return ---
if useGPU
    S_filtered = gather(S_filtered);
end

IQ_filtered = reshape(S_filtered, nZ, nX, nFrames);

fprintf('    SVD (rsvd k=%d): removed SV 1:%d of %d (%.2fs)\n', ...
    k, cutoffLow, nFrames, toc(tSVD));

end
