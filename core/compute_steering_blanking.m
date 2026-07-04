function info = compute_steering_blanking(steerAngle_deg, nRx, pitch_mm, c, fs_MHz, margin, voltage_pct)
% COMPUTE_STEERING_BLANKING  Compute RF sample blanking for steered plane waves.
%
% When a plane wave is steered at angle theta, elements fire sequentially
% with a delay spread of (nRx-1)*pitch*|sin(theta)|/c. During this window,
% electrical cross-talk between adjacent TX elements corrupts the earliest
% RF samples asymmetrically for +/- PI pulses, degrading PI cancellation.
%
% At high TX voltage (75-90%), amplifier nonlinearity extends the
% contamination beyond the theoretical firing window. The margin is
% automatically increased for high-power acquisitions.
%
% Inputs:
%   steerAngle_deg - Steering angle [degrees]
%   nRx            - Number of receive elements
%   pitch_mm       - Element pitch [mm]
%   c              - Speed of sound [m/s]
%   fs_MHz         - Sampling frequency [MHz]
%   margin         - Base safety margin (default 1.5). Set to [] for auto.
%   voltage_pct    - TX voltage percentage (default []). If provided,
%                    margin is automatically scaled:
%                      0-25%:  margin * 1.0
%                      25-50%: margin * 1.2
%                      50-75%: margin * 1.5
%                      75-100%: margin * 2.0
%
% Output:
%   info - Struct with fields:
%     .nBlank             - Samples to blank
%     .minDepth_mm        - Minimum valid imaging depth
%     .delaySpread_us     - TX delay spread
%     .delaySpread_samples - Delay spread in samples
%     .angle_deg          - Input angle
%     .margin             - Final margin used (after voltage adjustment)
%     .voltage_pct        - TX voltage used

if nargin < 6 || isempty(margin), margin = 1.5; end
if nargin < 7, voltage_pct = []; end

% Adjust margin based on TX voltage
if ~isempty(voltage_pct)
    if voltage_pct > 75
        voltageScale = 2.0;
    elseif voltage_pct > 50
        voltageScale = 1.5;
    elseif voltage_pct > 25
        voltageScale = 1.2;
    else
        voltageScale = 1.0;
    end
    margin = margin * voltageScale;
end

c_mm_us = c * 1e-3;
theta_rad = abs(steerAngle_deg) * pi / 180;

% TX delay spread across the aperture
delaySpread_us = (nRx - 1) * pitch_mm * sin(theta_rad) / c_mm_us;

% Samples to blank (with margin)
nBlank = ceil(delaySpread_us * fs_MHz * margin);

% Minimum valid one-way depth
minDepth_mm = c_mm_us * delaySpread_us * margin / 2;

info.nBlank = nBlank;
info.minDepth_mm = minDepth_mm;
info.delaySpread_us = delaySpread_us;
info.delaySpread_samples = delaySpread_us * fs_MHz;
info.angle_deg = steerAngle_deg;
info.margin = margin;
info.voltage_pct = voltage_pct;

end
