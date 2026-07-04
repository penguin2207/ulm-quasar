function [x_refit, info] = quasar_refit(y, psf, x_fista, opts)
% QUASAR_REFIT  Post-LASSO (relaxed-LASSO at phi=0) least-squares amplitude
% refit on the FISTA support.
%
% The core of QUASAR: L1 soft-thresholding in FISTA compresses the
% amplitude dynamic range because the threshold shrinks all values
% toward zero. The post-LASSO refit removes that shrinkage:
%   1. Use FISTA solution to detect support S = {i : x_fista(i) > 0}
%   2. Solve unregularized LS restricted to S:
%      x_S = argmin ||y - A_S * x_S||^2
%      where A_S is the PSF operator restricted to support pixels
%
% This recovers ~unbiased amplitude estimates that preserve the
% concentration-dependent scaling needed for quantitative imaging.
%
% NAMING (per the 2026-06-08 correction): this is POST-LASSO, i.e. relaxed
% LASSO at phi=0 (Belloni & Chernozhukov 2013; Meinshausen 2007): L1 selects
% the support, then unregularized LS re-estimates the amplitudes on it. It is
% NOT 'debiased'/'desparsified' LASSO (Zhang & Zhang 2014; van de Geer et al.
% 2014; Javanmard & Montanari 2014), which is a one-step correction of the
% FULL coefficient vector for inference (a different procedure). Use
% post-LASSO/relaxed-LASSO in the writeup, not debiased LASSO.
%
% Inputs:
%   y        - [nZ x nX] observed power image
%   psf      - [nZ x nX] PSF (same size, centered)
%   x_fista  - [nZ x nX] FISTA sparse solution (used for support detection)
%   opts     - Struct with optional parameters:
%     .supportThresh - Fraction of max to threshold support (default 0)
%                      Set >0 to prune weak detections before refit
%     .maxIterCG     - Max conjugate gradient iterations (default 50)
%     .tolCG         - CG convergence tolerance (default 1e-6)
%     .useGPU        - Use GPU (default false)
%     .verbose       - Print progress (default false)
%
% Outputs:
%   x_refit  - [nZ x nX] refit (post-LASSO) amplitude image
%   info     - Struct with:
%     .supportSize  - Number of support pixels
%     .residualNorm - ||y - A*x_refit|| after refit
%     .amplitudeRatio - mean(x_refit(S)) / mean(x_fista(S)) (bias measure)
%     .time         - Computation time

if nargin < 4, opts = struct(); end

supportThresh = getfield_default(opts, 'supportThresh', 0);
maxIterCG     = getfield_default(opts, 'maxIterCG', 50);
tolCG         = getfield_default(opts, 'tolCG', 1e-6);
useGPU        = getfield_default(opts, 'useGPU', false);
verbose       = getfield_default(opts, 'verbose', false);

tStart = tic;

[nZ, nX] = size(y);

% Step 1: Detect support from FISTA solution
if supportThresh > 0
    thresh = supportThresh * max(x_fista(:));
    support = x_fista > thresh;
else
    support = x_fista > 0;
end
supportIdx = find(support);
nSupport = numel(supportIdx);

if verbose
    fprintf('  QUASAR refit: %d support pixels (%.2f%% of grid)\n', ...
        nSupport, 100*nSupport/(nZ*nX));
end

if nSupport == 0
    x_refit = zeros(nZ, nX, 'single');
    info.supportSize = 0;
    info.residualNorm = norm(y(:));
    info.amplitudeRatio = NaN;
    info.time = toc(tStart);
    return;
end

% Step 2: Solve LS restricted to support using conjugate gradient
% We solve: A_S' * A_S * x_S = A_S' * y
% where A_S applies PSF convolution then masks to support

if useGPU
    y = gpuArray(single(y));
    H = gpuArray(single(fft2(ifftshift(psf))));
else
    y = single(y);
    H = single(fft2(ifftshift(psf)));
end
Hconj = conj(H);

% Operators restricted to support
% A_S(x_S): zero-pad to full grid, convolve with PSF
% A_S'(r): convolve with PSF', mask to support
applyA_full = @(x) real(ifft2(H .* fft2(x)));
applyAT_full = @(x) real(ifft2(Hconj .* fft2(x)));

% Mask operators (support <-> full grid)
fullToSupport = @(x) x(supportIdx);
supportToFull = @(xs) embedSupport(xs, supportIdx, nZ, nX, useGPU);

% Normal equation operator: A_S' * A_S
applyNormal = @(xs) fullToSupport(applyAT_full(applyA_full(supportToFull(xs))));

% Right-hand side: A_S' * y
rhs = fullToSupport(applyAT_full(y));

% Conjugate gradient solve -- warm-start from FISTA solution on support
% (instead of zeros, so CG only needs to correct the L1 bias)
if useGPU
    x_s = gpuArray(single(x_fista(supportIdx)));
else
    x_s = single(x_fista(supportIdx));
end
r = rhs - applyNormal(x_s);
p = r;
rsold = r' * r;

for iter = 1:maxIterCG
    Ap = applyNormal(p);
    alpha = rsold / (p' * Ap + eps);
    x_s = x_s + alpha * p;
    r = r - alpha * Ap;
    rsnew = r' * r;
    
    if sqrt(rsnew) < tolCG * sqrt(rhs' * rhs + eps)
        if verbose
            fprintf('  CG converged at iter %d\n', iter);
        end
        break;
    end
    
    p = r + (rsnew / (rsold + eps)) * p;
    rsold = rsnew;
end

% Enforce non-negativity (physical constraint: power is non-negative)
x_s = max(x_s, 0);

% Embed back into full grid
x_refit = supportToFull(x_s);

if useGPU
    x_refit = gather(x_refit);
    x_fista_g = x_fista;
else
    x_fista_g = x_fista;
end

% Compute stats
residual = applyA_full(supportToFull(x_s)) - y;
resNorm = gather(norm(residual(:)));

% Amplitude ratio: how much did the refit change amplitudes?
fista_on_support = x_fista_g(supportIdx);
refit_on_support = gather(x_s);
amplitudeRatio = mean(refit_on_support) / (mean(fista_on_support(:)) + eps);

info.supportSize = nSupport;
info.residualNorm = resNorm;
info.amplitudeRatio = amplitudeRatio;
info.cgIters = min(iter, maxIterCG);
info.time = toc(tStart);

if verbose
    fprintf('  Refit: support=%d, residual=%.2f, amp ratio=%.2f (>1 = FISTA was biased low)\n', ...
        nSupport, resNorm, amplitudeRatio);
end

end


function xfull = embedSupport(xs, idx, nZ, nX, useGPU)
    if useGPU
        xfull = gpuArray(zeros(nZ, nX, 'single'));
    else
        xfull = zeros(nZ, nX, 'single');
    end
    xfull(idx) = xs;
end

function val = getfield_default(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end
