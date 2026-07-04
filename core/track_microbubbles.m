function tracks = track_microbubbles(localizations, params, timestamps)
% TRACK_MICROBUBBLES  Frame-to-frame tracking via Kalman filter + Hungarian.
%
% Links microbubble localizations across frames into continuous tracks
% using a constant-velocity Kalman filter for state prediction and the
% Hungarian algorithm (Kuhn-Munkres) for optimal assignment.
%
% References:
%   Tang S, Song P, Trzasko JD, et al. "Kalman filter-based microbubble
%   tracking for robust super-resolution ultrasound microvessel imaging."
%   IEEE Trans UFFC. 2020;67(9):1738-1751.
%
%   Hingot V et al. used Hungarian method via simpletracker.
%   Sci Rep. 2019;9:2456.
%
% Inputs:
%   localizations - [N x 4] array: [x_mm, z_mm, amplitude, frameIdx]
%   params        - Struct with tracking parameters:
%                   .maxDisp_mm     - Max displacement per frame [mm]
%                   .maxGapFrames   - Max frames for gap closing
%                   .minTrackLength - Minimum track length to retain
%                   .kalman.processNoise  - Process noise variance
%                   .kalman.measNoise     - Measurement noise variance [mm]
%   timestamps    - [nFrames x 1] timestamps [ms] (for velocity estimation)
%
% Output:
%   tracks - Cell array. Each cell is [M x 4]: [x_mm, z_mm, amplitude, frameIdx]

if isempty(localizations)
    tracks = {};
    return;
end

% Get frame indices
frames = unique(localizations(:, 4));
nFrames = numel(frames);

fprintf('    Tracking: %d localizations across %d frames\n', ...
    size(localizations, 1), nFrames);

% --- Initialize data structures ---
% Active track list: each track has a Kalman filter state and history
activeTracksState = {};   % Kalman states: {[x, z, vx, vz]}
activeTracksP     = {};   % Covariance matrices
activeTracksHist  = {};   % History: [x, z, amp, frame]
activeTracksAge   = [];   % Frames since last assignment
nextTrackId = 1;

% Completed tracks
completedTracks = {};

% Kalman filter matrices (constant velocity model in 2D)
dt = 1;  % Normalized time step (1 frame)
F = [1 0 dt 0;    % State transition
     0 1 0  dt;
     0 0 1  0;
     0 0 0  1];
H = [1 0 0 0;     % Measurement matrix (observe position only)
     0 1 0 0];

q = params.kalman.processNoise;
Q = q * [dt^3/3  0      dt^2/2  0;      % Process noise
         0       dt^3/3  0      dt^2/2;
         dt^2/2  0       dt      0;
         0       dt^2/2  0       dt];
     
r = params.kalman.measNoise^2;
R = r * eye(2);  % Measurement noise

% --- Frame-by-frame processing ---
% --- Frame-by-frame processing ---
wb = [];
if nFrames > 100
    wb = waitbar(0, 'Tracking microbubbles...', 'Name', 'LAT-ULM Tracking');
end

for iFrame = 1:nFrames
    if ~isempty(wb) && mod(iFrame, 200) == 0 && ishandle(wb)
        waitbar(iFrame/nFrames, wb, sprintf(...
            'Frame %d/%d | %d active tracks | %d completed', ...
            iFrame, nFrames, numel(activeTracksState), numel(completedTracks)));
    end
    
    frameIdx = frames(iFrame);
    
    % Get localizations in this frame
    mask = localizations(:, 4) == frameIdx;
    frameLocs = localizations(mask, :);  % [nDets x 4]
    nDets = size(frameLocs, 1);
    nActive = numel(activeTracksState);
    
    if nActive == 0 && nDets == 0
        continue;
    end
    
    % --- Kalman prediction for all active tracks ---
    predictedPos = zeros(nActive, 2);
    for t = 1:nActive
        activeTracksState{t} = F * activeTracksState{t};
        activeTracksP{t}     = F * activeTracksP{t} * F' + Q;
        predictedPos(t, :)   = activeTracksState{t}(1:2)';
    end
    
    % --- Build cost matrix (distance between predictions and detections) ---
    if nActive > 0 && nDets > 0
        costMatrix = zeros(nActive, nDets);
        for t = 1:nActive
            for d = 1:nDets
                dist = sqrt((predictedPos(t,1) - frameLocs(d,1))^2 + ...
                           (predictedPos(t,2) - frameLocs(d,2))^2);
                if dist > params.maxDisp_mm
                    costMatrix(t, d) = Inf;  % Too far, disallow assignment
                else
                    costMatrix(t, d) = dist;
                end
            end
        end
        
        % --- Hungarian algorithm for optimal assignment ---
        [assignment, ~] = hungarian_algorithm(costMatrix);
        
    else
        assignment = zeros(max(nActive, 1), 1);
    end
    
    % --- Update assigned tracks ---
    assignedDets = false(nDets, 1);
    
    for t = 1:nActive
        if assignment(t) > 0 && ~isinf(costMatrix(t, assignment(t)))
            d = assignment(t);
            assignedDets(d) = true;
            
            % Kalman update
            z_meas = frameLocs(d, 1:2)';
            y = z_meas - H * activeTracksState{t};  % Innovation
            S = H * activeTracksP{t} * H' + R;
            K = activeTracksP{t} * H' / S;  % Kalman gain
            
            activeTracksState{t} = activeTracksState{t} + K * y;
            activeTracksP{t}     = (eye(4) - K * H) * activeTracksP{t};
            activeTracksAge(t)   = 0;
            
            % Append to history
            activeTracksHist{t} = [activeTracksHist{t}; ...
                activeTracksState{t}(1), activeTracksState{t}(2), ...
                frameLocs(d, 3), frameIdx];
        else
            % Track not assigned this frame
            activeTracksAge(t) = activeTracksAge(t) + 1;
        end
    end
    
    % --- Terminate old tracks ---
    terminate = activeTracksAge > params.maxGapFrames;
    for t = find(terminate)
        hist = activeTracksHist{t};
        if ~isempty(hist) && isnumeric(hist) && size(hist, 1) >= params.minTrackLength
            completedTracks{end+1} = hist; %#ok<AGROW>
        end
    end
    activeTracksState(terminate) = [];
    activeTracksP(terminate)     = [];
    activeTracksHist(terminate)  = [];
    activeTracksAge(terminate)   = [];
    
    % --- Start new tracks for unassigned detections ---
    for d = find(~assignedDets)'
        newState = [frameLocs(d, 1); frameLocs(d, 2); 0; 0];  % [x, z, vx, vz]
        newP = blkdiag(R, 10 * r * eye(2));  % Large initial velocity uncertainty
        
        activeTracksState{end+1} = newState; %#ok<AGROW>
        activeTracksP{end+1}     = newP; %#ok<AGROW>
        activeTracksHist{end+1}  = [frameLocs(d, 1:3), frameIdx]; %#ok<AGROW>
        activeTracksAge(end+1)   = 0; %#ok<AGROW>
    end
end

% --- Finalize remaining active tracks ---
for t = 1:numel(activeTracksState)
    hist = activeTracksHist{t};
    if ~isempty(hist) && isnumeric(hist) && size(hist, 1) >= params.minTrackLength
        completedTracks{end+1} = hist; %#ok<AGROW>
    end
end

tracks = completedTracks;

if ~isempty(wb) && ishandle(wb), close(wb); end

fprintf('    Completed: %d tracks (min length = %d)\n', numel(tracks), params.minTrackLength);

end
