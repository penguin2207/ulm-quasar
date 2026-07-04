function [subRow, subCol] = intensity_weighted_centroid(roi)
% INTENSITY_WEIGHTED_CENTROID  Sub-pixel localization via intensity weighting.
%
% Computes the intensity-weighted centroid of a region of interest to
% achieve sub-pixel localization precision. This is the simplest and most
% widely used localization method in ULM.
%
% Reference: Desailly Y, Couture O, Fink M, Tanter M. "Sono-activated
%            ultrasound localization microscopy."
%            Appl Phys Lett. 2013;103(17):174107.
%
% Also: Heiles et al. (PALA) benchmark shows weighted centroid achieves
%       comparable precision to Gaussian fitting for high-SNR bubbles.
%       Nat Biomed Eng. 2022;6(5):605-616.
%
% Input:
%   roi - [nR x nC] intensity image patch around a detected bubble
%
% Outputs:
%   subRow - Sub-pixel row position within ROI (1-indexed)
%   subCol - Sub-pixel column position within ROI (1-indexed)

[nR, nC] = size(roi);

% Background subtraction: remove minimum to avoid bias from background
roi = roi - min(roi(:));

% Ensure non-negative (should already be from envelope data)
roi = max(roi, 0);

totalIntensity = sum(roi(:));

if totalIntensity == 0
    % Fallback to center of ROI
    subRow = (nR + 1) / 2;
    subCol = (nC + 1) / 2;
    return;
end

% Intensity-weighted centroid
[colGrid, rowGrid] = meshgrid(1:nC, 1:nR);
subRow = sum(rowGrid(:) .* roi(:)) / totalIntensity;
subCol = sum(colGrid(:) .* roi(:)) / totalIntensity;

end
