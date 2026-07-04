function [psf, psfParams] = build_sushi_psf(freq_MHz, c, pitch_mm, nRx, ...
    srPixel_mm, srGridSize, steerAngle_deg)
% BUILD_SUSHI_PSF  Construct analytical PSF model for SUSHI sparse recovery.
%
% Builds a 2D point spread function based on probe parameters for use
% in the SUSHI deconvolution formulation.
%
% The PSF is modeled as:
%   Axial:   sinc envelope from pulse bandwidth
%   Lateral: sinc from receive aperture (plane wave, no TX focus)
%
% Inputs:
%   freq_MHz       - TX center frequency [MHz]
%   c              - Speed of sound [m/s]
%   pitch_mm       - Element pitch [mm]
%   nRx            - Number of receive elements
%   srPixel_mm     - Super-resolution pixel size [mm]
%   srGridSize     - [nZ, nX] of the SR grid
%   steerAngle_deg - Steering angle [degrees] (0 for standard)
%
% Outputs:
%   psf        - [nZ x nX] normalized PSF (same size as SR grid, centered)
%   psfParams  - Struct with PSF parameters for documentation

lambda_mm = c * 1e-3 / freq_MHz;  % Wavelength [mm]
aperture_mm = (nRx - 1) * pitch_mm;  % Receive aperture [mm]

% PSF dimensions
% Axial FWHM ~ lambda (for ~1 cycle pulse)
% Lateral FWHM ~ lambda * F# (F# = depth/aperture, varies with depth)
% Use a representative depth (center of typical FOV)
typicalDepth_mm = 20;  % Representative depth (tubes at ~18-22 mm in phantom)
Fnum = typicalDepth_mm / aperture_mm;

axialFWHM_mm = lambda_mm;           % ~1 wavelength for broadband pulse
lateralFWHM_mm = lambda_mm * Fnum;  % Diffraction-limited lateral

% Convert to pixels
axialFWHM_px = axialFWHM_mm / srPixel_mm;
lateralFWHM_px = lateralFWHM_mm / srPixel_mm;

% Build 2D PSF (Gaussian approximation of sinc mainlobe)
% sigma = FWHM / (2*sqrt(2*ln(2)))
sigmaZ = axialFWHM_px / 2.355;
sigmaX = lateralFWHM_px / 2.355;

% Create PSF kernel (odd-sized, centered)
kernelHalfZ = ceil(3 * sigmaZ);
kernelHalfX = ceil(3 * sigmaX);
[kx, kz] = meshgrid(-kernelHalfX:kernelHalfX, -kernelHalfZ:kernelHalfZ);
psfKernel = exp(-0.5 * (kz.^2/sigmaZ^2 + kx.^2/sigmaX^2));
psfKernel = psfKernel / sum(psfKernel(:));  % Normalize to unit sum

% Embed in full SR grid (zero-padded, centered)
psf = zeros(srGridSize);
centerZ = ceil(srGridSize(1)/2);
centerX = ceil(srGridSize(2)/2);
z1 = max(1, centerZ - kernelHalfZ);
z2 = min(srGridSize(1), centerZ + kernelHalfZ);
x1 = max(1, centerX - kernelHalfX);
x2 = min(srGridSize(2), centerX + kernelHalfX);

kz1 = kernelHalfZ + 1 - (centerZ - z1);
kz2 = kernelHalfZ + 1 + (z2 - centerZ);
kx1 = kernelHalfX + 1 - (centerX - x1);
kx2 = kernelHalfX + 1 + (x2 - centerX);

psf(z1:z2, x1:x2) = psfKernel(kz1:kz2, kx1:kx2);

% Store parameters
psfParams.lambda_mm = lambda_mm;
psfParams.aperture_mm = aperture_mm;
psfParams.axialFWHM_mm = axialFWHM_mm;
psfParams.lateralFWHM_mm = lateralFWHM_mm;
psfParams.axialFWHM_px = axialFWHM_px;
psfParams.lateralFWHM_px = lateralFWHM_px;
psfParams.kernelSize = size(psfKernel);
psfParams.typicalDepth_mm = typicalDepth_mm;
psfParams.Fnum = Fnum;

end
