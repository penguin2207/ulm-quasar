function tracksQC = filter_tracks_quality(tracks, trackParams)
% FILTER_TRACKS_QUALITY  Kinematic-coherence track gate (default pass-through).
%
%   tracksQC = filter_tracks_quality(tracks, trackParams)
%
% Keeps only tracks passing the QC thresholds in trackParams, when set:
%   .minMeanAmp      - drop tracks with mean amplitude below this
%   .minStraightness - drop tracks with net/path displacement below this
%                      (flow-coherence; clutter coincidences are low, ~<0.5)
%   .maxGapInTrack   - drop tracks whose largest frame gap exceeds this
% When NONE are set (PROFILE_APR17 / reproduction), returns tracks unchanged.
%
% This is the lever (not amplitude threshold) that separates sparse low-conc
% bubbles from clutter; see FINDINGS_2026_06_11_lowconc_floor.md.
%
% NOTE: assumes each track is a numeric [nPts x >=4] array with columns
% [x_mm, z_mm, amp, frameIdx] (the localization rows track_microbubbles links).
% VERIFY this column order against track_microbubbles output before relying on
% the straightness/amp gates; the pass-through path (apr17) does not touch it.

    tracksQC = tracks(:);
    if isempty(tracksQC), return; end

    hasAmp = isfield(trackParams,'minMeanAmp')      && ~isempty(trackParams.minMeanAmp);
    hasStr = isfield(trackParams,'minStraightness') && ~isempty(trackParams.minStraightness);
    hasGap = isfield(trackParams,'maxGapInTrack')   && ~isempty(trackParams.maxGapInTrack);
    if ~(hasAmp || hasStr || hasGap)
        return;   % pass-through (no QC requested)
    end

    keep = true(numel(tracksQC), 1);
    for i = 1:numel(tracksQC)
        T = tracksQC{i};
        if ~isnumeric(T) || size(T,1) < 2, continue; end
        P = T(:, 1:2);
        if hasAmp && size(T,2) >= 3 && mean(T(:,3)) < trackParams.minMeanAmp
            keep(i) = false; continue;
        end
        if hasStr
            pathLen = sum(sqrt(sum(diff(P).^2, 2)));
            straight = norm(P(end,:) - P(1,:)) / max(pathLen, eps);
            if straight < trackParams.minStraightness
                keep(i) = false; continue;
            end
        end
        if hasGap && size(T,2) >= 4
            if max(diff(T(:,4))) > trackParams.maxGapInTrack
                keep(i) = false; continue;
            end
        end
    end
    tracksQC = tracksQC(keep);
end
