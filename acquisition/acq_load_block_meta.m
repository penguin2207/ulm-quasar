function meta = acq_load_block_meta(folder, basename, modeName, sosOverride, pitchOverride)
% ACQ_LOAD_BLOCK_META  Load probe and acquisition metadata from one VADA block.
%
% Reads the first eventsPerFrame events of a block to determine:
%   - probe parameters (elements, pitch, fs)
%   - acquisition pattern (events/frame, angles, polarities)
%   - angle/polarity pairing (anglePairs struct)
%   - frame rate, depth offset, sample count
%
% Inputs:
%   folder        path to .vada folder
%   basename      base filename (no extension)
%   modeName      file extension ('.vada')
%   sosOverride   speed of sound override [m/s], default 1540
%   pitchOverride pitch override [mm] if metadata reports 0, default 0.300
%
% Output: meta struct with fields:
%   .filename
%   .probeName, .nElements, .pitch_mm
%   .fs_MHz, .depthOffset_mm, .nSamples
%   .c                          (speed of sound, m/s)
%   .txFreq_MHz, .lambda_mm
%   .eventsPerFrame, .numAngles, .hasPI
%   .anglePairs                 struct array per unique angle
%   .zeroAngleIdx, .zeroEventIdx (set if 0-deg angle present)
%   .totalEvents, .numCompoundFrames
%   .frameRate_Hz
%   .elemPos_mm
%   .txVoltage                  high voltage rail percent
%
% This is a "fast" metadata read: ~30 events of one block, no full data load.

if nargin < 3, modeName      = '.vada'; end
if nargin < 4, sosOverride   = [];      end  % empty = use XML value
if nargin < 5, pitchOverride = 0.300;   end

meta = struct();
meta.filename = basename;

% --- Probe events to detect events/frame ---
numProbe = 30;
[VadaProbe, Param, TxrParam, BlockConfig] = VsiVadaDataRead( ...
    folder, basename, 1:numProbe, modeName);

% Parse angles + polarities for each probe event
probeAngles = zeros(numProbe, 1);
probePolar  = zeros(numProbe, 1);
for ev = 1:numProbe
    if isfield(VadaProbe(ev).TxDelay, 'angle')
        probeAngles(ev) = VadaProbe(ev).TxDelay.angle;
    end
    if isfield(VadaProbe(ev).Waveform, 'Channel') && ...
            isfield(VadaProbe(ev).Waveform.Channel(1), 'invert')
        probePolar(ev) = VadaProbe(ev).Waveform.Channel(1).invert;
    end
end

% Find compound-frame period via repeating signature
probeSig = probeAngles * 10 + probePolar;
eventsPerFrame = 0;

% Short-circuit for truly single-event sequences (e.g. L38xp April 8-9,
% which is 1 plane wave + no PI + no compounding, so every event has
% identical signature). The pattern-detection loop below cannot detect
% patLen=1 because it skips ev=2 to avoid a degenerate match, so this
% branch catches the all-same-signature case explicitly.
if numel(unique(probeSig(1:numProbe))) == 1
    eventsPerFrame = 1;
else
    for ev = 2:numProbe
        if probeSig(ev) == probeSig(1) && ev > 2
            patLen = ev - 1;
            if ev + patLen - 1 <= numProbe
                if all(probeSig(ev:ev+patLen-1) == probeSig(1:patLen))
                    eventsPerFrame = patLen;
                    break;
                end
            end
        end
    end
end

if eventsPerFrame == 0
    % Fallback: pattern detection failed; assume single event
    eventsPerFrame = 1;
end

% --- Detect angles and PI ---
firstFrameAngles = probeAngles(1:eventsPerFrame);
firstFramePolar  = probePolar(1:eventsPerFrame);
uniqueAngles = unique(firstFrameAngles, 'stable');
hasPI = true;
for a = 1:numel(uniqueAngles)
    if sum(firstFrameAngles == uniqueAngles(a)) < 2
        hasPI = false;
        break;
    end
end

% Build angle pairs
anglePairs = struct('angle', {}, 'posIdx', {}, 'negIdx', {}, 'rxElements', {});
for a = 1:numel(uniqueAngles)
    anglePairs(a).angle = uniqueAngles(a);
    idxs = find(firstFrameAngles == uniqueAngles(a));
    if numel(idxs) >= 2
        if firstFramePolar(idxs(1)) == 0
            anglePairs(a).posIdx = idxs(1); anglePairs(a).negIdx = idxs(2);
        else
            anglePairs(a).posIdx = idxs(2); anglePairs(a).negIdx = idxs(1);
        end
    else
        anglePairs(a).posIdx = idxs(1);
        anglePairs(a).negIdx = [];
    end
    anglePairs(a).rxElements = VadaProbe(anglePairs(a).posIdx).Elements;
end

% --- Probe geometry ---
rawPitch = TxrParam.ArrayPitch;
if rawPitch == 0 || isnan(rawPitch)
    pitch_mm = pitchOverride;
elseif rawPitch < 10
    pitch_mm = rawPitch;
else
    pitch_mm = rawPitch / 1000;
end

% --- Speed of sound override ---
c_xml = Param.SoSMedia;
if ~isempty(sosOverride)
    c = sosOverride;
elseif c_xml > 0
    c = c_xml;
else
    c = 1540;
end

% --- TX frequency ---
txFreq_MHz = 6;  % default
if isfield(VadaProbe(1).Waveform, 'Channel') && ...
        ~isempty(VadaProbe(1).Waveform.Channel) && ...
        isfield(VadaProbe(1).Waveform.Channel(1), 'frequency')
    txFreq_MHz = VadaProbe(1).Waveform.Channel(1).frequency;
end
lambda_mm = c * 1e-3 / txFreq_MHz;

% --- Sample count from first event ---
nSamples = size(VadaProbe(1).Data, 1);

% --- TX voltage from XML for blanking margin ---
xmlPath = fullfile(folder, [basename modeName '.xml']);
txVoltage = [];
if exist(xmlPath, 'file')
    try
        xmlParams = read_vada_xml_params(xmlPath);
        voltHi = get_param(xmlParams, 'Vada-Mode/Voltage-Rail-High', '');
        if ~isempty(voltHi)
            txVoltage = str2double(voltHi);
        end
    catch
        % silent: blanking will use base margin
    end
end

% --- Total events and compound frame count ---
totalEvents = numel(BlockConfig.PulseSequences(1).Events);
numCompoundFrames = floor(totalEvents / eventsPerFrame);

% --- Frame rate from PRI x events/frame ---
% Default PRI = 150 us (verified from XML for both probes Apr 8-9)
% Use actual PRI from XML if available; fall back to 150 us.
priUs = 150;  % microseconds, default
cfgPath = fullfile(folder, [basename modeName '.cfg.xml']);
if exist(cfgPath, 'file')
    try
        priFromXml = parse_pri_from_cfg(cfgPath);
        if ~isnan(priFromXml) && priFromXml > 0
            priUs = priFromXml;
        end
    catch
        % silent: use default
    end
end
frameRate_Hz = 1e6 / (eventsPerFrame * priUs);

% --- Element positions (centered on array) ---
elemPos_mm = ((1:TxrParam.ArrayNumElements) - (TxrParam.ArrayNumElements+1)/2) * pitch_mm;

% --- 0-degree angle index for zeroOnly mode ---
zeroAngleIdx = find([anglePairs.angle] == 0);
if isempty(zeroAngleIdx)
    [~, zeroAngleIdx] = min(abs([anglePairs.angle]));
end
zeroEventIdx = anglePairs(zeroAngleIdx).posIdx;

% --- Pack ---
meta.probeName        = TxrParam.Name;
meta.nElements        = TxrParam.ArrayNumElements;
meta.pitch_mm         = pitch_mm;
meta.fs_MHz           = Param.SampleFreq;
meta.depthOffset_mm   = Param.DepthOffset;
meta.nSamples         = nSamples;
meta.c                = c;
meta.txFreq_MHz       = txFreq_MHz;
meta.lambda_mm        = lambda_mm;
meta.eventsPerFrame   = eventsPerFrame;
meta.numAngles        = numel(uniqueAngles);
meta.hasPI            = hasPI;
meta.anglePairs       = anglePairs;
meta.zeroAngleIdx     = zeroAngleIdx;
meta.zeroEventIdx     = zeroEventIdx;
meta.totalEvents      = totalEvents;
meta.numCompoundFrames = numCompoundFrames;
meta.frameRate_Hz     = frameRate_Hz;
meta.elemPos_mm       = elemPos_mm;
meta.txVoltage        = txVoltage;
meta.priUs            = priUs;

end


function pri_us = parse_pri_from_cfg(cfgPath)
% Quick regex extraction of PRI from .vada.cfg.xml
pri_us = NaN;
try
    txt = fileread(cfgPath);
    tok = regexp(txt, 'pri="([\d.]+)"', 'tokens', 'once');
    if ~isempty(tok)
        pri_us = str2double(tok{1});
    end
catch
    % silent
end
end
