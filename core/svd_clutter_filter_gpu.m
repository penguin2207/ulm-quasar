function IQ_filtered = svd_clutter_filter_gpu(IQ_stack, cutoffLow, cutoffHigh, useGPU)
% SVD_CLUTTER_FILTER_GPU  GPU-accelerated spatiotemporal SVD clutter rejection.
%
% Decomposes spatiotemporal IQ stack via SVD, removes tissue clutter
% (high singular values) and optionally noise (low singular values).
%
% Reference: Demene C, Deffieux T, Pernot M, et al.
%            IEEE Trans Med Imaging. 2015;34(11):2271-2285.
%
% Inputs:
%   IQ_stack   - [nZ x nX x nFrames] complex IQ image stack
%   cutoffLow  - Remove singular vectors 1:cutoffLow (tissue)
%   cutoffHigh - Remove above cutoffHigh (noise). [] = keep all
%   useGPU     - Boolean: use GPU acceleration
%
% Output:
%   IQ_filtered - [nZ x nX x nFrames] filtered stack

[nZ, nX, nFrames] = size(IQ_stack);
nPixels = nZ * nX;

tSVD = tic;

% Reshape to Casorati matrix [pixels x frames]
S = reshape(IQ_stack, nPixels, nFrames);

% Move to GPU if requested
if useGPU
    S = gpuArray(single(S));
end

% SVD (MATLAB's svd works on gpuArray automatically)
% For large matrices, economy SVD is sufficient
if nPixels > nFrames
    % Typical case: more pixels than frames
    [U, Sigma, V] = svd(S, 'econ');
else
    [U, Sigma, V] = svd(S, 'econ');
end

% Determine keep range
svStart = cutoffLow + 1;
if isempty(cutoffHigh) || cutoffHigh == 0
    svEnd = min(size(Sigma));
else
    svEnd = min(cutoffHigh, min(size(Sigma)));
end

% Reconstruct with selected components
S_filtered = U(:, svStart:svEnd) * Sigma(svStart:svEnd, svStart:svEnd) * V(:, svStart:svEnd)';

% Gather from GPU and reshape
if useGPU
    S_filtered = gather(S_filtered);
end

IQ_filtered = reshape(S_filtered, nZ, nX, nFrames);

fprintf('    SVD: removed SV 1:%d, kept %d:%d of %d (%.1fs)\n', ...
    cutoffLow, svStart, svEnd, min(size(Sigma)), toc(tSVD));

end
