function [x, info] = sushi_sparse_recovery(y, psf, lambda, method, opts)
% SUSHI_SPARSE_RECOVERY  Sparse recovery for SUSHI super-resolution.
%
% Solves the deconvolution problem: minimize (1/2)||y - H*x||^2 + lambda*||x||_1
% where H is the PSF convolution operator (applied via FFT).
%
% Methods:
%   'fista'     - Fast Iterative Shrinkage-Thresholding (standard SUSHI)
%   'ista'      - Basic ISTA (for comparison, slower convergence)
%
% Inputs:
%   y        - [nZ x nX] observed power image (real, non-negative)
%   psf      - [nZ x nX] PSF (same size as y, centered, from build_sushi_psf)
%   lambda   - L1 regularization weight (higher = sparser, less noisy)
%   method   - 'fista' (default) or 'ista'
%   opts     - Struct with optional parameters:
%     .maxIter  - Maximum iterations (default 100)
%     .tol      - Convergence tolerance on relative change (default 1e-4)
%     .verbose  - Print progress (default false)
%     .nonNeg   - Enforce non-negativity (default true, physical for power)
%     .useGPU   - Use GPU arrays (default false)
%
% Outputs:
%   x        - [nZ x nX] sparse super-resolution image
%   info     - Struct with convergence info:
%     .cost     - Cost function at each iteration
%     .nIter    - Iterations used
%     .support  - Number of nonzero pixels
%     .maxVal   - Maximum recovered value
%     .time     - Computation time

if nargin < 4 || isempty(method), method = 'fista'; end
if nargin < 5, opts = struct(); end

% Default options
maxIter = getfield_default(opts, 'maxIter', 100);
tol     = getfield_default(opts, 'tol', 1e-4);
verbose = getfield_default(opts, 'verbose', false);
nonNeg  = getfield_default(opts, 'nonNeg', true);
useGPU  = getfield_default(opts, 'useGPU', false);

tStart = tic;

% Precompute PSF in frequency domain for fast convolution
if useGPU
    y = gpuArray(single(y));
    H = gpuArray(single(fft2(ifftshift(psf))));  % Centered PSF -> FFT
else
    y = single(y);
    H = single(fft2(ifftshift(psf)));
end
Hconj = conj(H);

% Forward model: A(x) = real(ifft2(H .* fft2(x)))
applyA   = @(x) real(ifft2(H .* fft2(x)));
applyAT  = @(x) real(ifft2(Hconj .* fft2(x)));

% Lipschitz constant: L = max(|H|^2)
% Step size = 1/L
L = max(abs(H(:)).^2);
stepSize = 1 / L;

% Compute lambda_max: the smallest lambda that gives x=0.
% At x=0, gradient = A^T(A*0 - y) = -A^T*y, so lambda_max = max(|A^T*y|) / stepSize^{-1}
% = max(|A^T*y|) * stepSize.  Any lambda >= lambda_max produces all-zero solution.
% User-supplied lambda is treated as a FRACTION of lambda_max (0 < lambda < 1).
ATy = applyAT(y);
if nonNeg
    lambdaMax = max(ATy(:)) * stepSize;  % Only positive gradient matters
else
    lambdaMax = max(abs(ATy(:))) * stepSize;
end
lambdaAbs = lambda * lambdaMax;

% Soft-thresholding operator
thresh = lambdaAbs * stepSize;
if nonNeg
    softThresh = @(z) max(z - thresh, 0);  % Non-negative soft threshold
else
    softThresh = @(z) sign(z) .* max(abs(z) - thresh, 0);
end

% Initialize
[nZ, nX] = size(y);
x = zeros(nZ, nX, 'like', y);  % Start from zero
z = x;  % FISTA auxiliary variable
t = 1;  % FISTA momentum parameter

costHistory = zeros(maxIter, 1);

for iter = 1:maxIter
    % Gradient step
    residual = applyA(z) - y;
    grad = applyAT(residual);
    z_grad = z - stepSize * grad;
    
    % Proximal step (soft thresholding)
    x_new = softThresh(z_grad);
    
    % Cost function
    dataFit = 0.5 * sum(residual(:).^2);
    l1Norm = lambdaAbs * sum(abs(x_new(:)));
    costHistory(iter) = dataFit + l1Norm;
    
    if strcmp(method, 'fista')
        % FISTA momentum update
        t_new = (1 + sqrt(1 + 4*t^2)) / 2;
        z = x_new + ((t - 1) / t_new) * (x_new - x);
        t = t_new;
    else
        % ISTA: no momentum
        z = x_new;
    end
    
    % Convergence check
    if iter > 1
        relChange = abs(costHistory(iter) - costHistory(iter-1)) / (abs(costHistory(iter-1)) + eps);
        if relChange < tol
            if verbose
                fprintf('  Converged at iter %d (relChange=%.2e)\n', iter, relChange);
            end
            break;
        end
    end
    
    x = x_new;
    
    if verbose && mod(iter, 20) == 0
        fprintf('  Iter %d: cost=%.4e, support=%d, max=%.4f\n', ...
            iter, costHistory(iter), nnz(x), max(x(:)));
    end
end

% Gather from GPU if needed
if useGPU
    x = gather(x);
    costHistory = gather(costHistory);
end

% Output info
info.cost = costHistory(1:iter);
info.nIter = iter;
info.support = nnz(x);
info.maxVal = max(x(:));
info.time = toc(tStart);
info.lambda = lambda;
info.lambdaAbs = gather(lambdaAbs);
info.lambdaMax = gather(lambdaMax);
info.method = method;
info.stepSize = gather(stepSize);

end


function val = getfield_default(s, field, default)
    if isfield(s, field)
        val = s.(field);
    else
        val = default;
    end
end
