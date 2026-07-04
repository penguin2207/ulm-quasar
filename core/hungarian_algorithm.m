function [assignment, cost] = hungarian_algorithm(costMatrix)
% HUNGARIAN_ALGORITHM  Kuhn-Munkres optimal assignment.
%
% Solves the rectangular assignment problem: given a cost matrix where
% entry (i,j) is the cost of assigning track i to detection j, find the
% assignment that minimizes total cost.
%
% Handles Inf entries (forbidden assignments) and rectangular matrices
% (more tracks than detections or vice versa).
%
% Used in ULM for frame-to-frame microbubble association.
% Reference: Kuhn HW (1955). "The Hungarian method for the assignment
%            problem." Naval Research Logistics Quarterly, 2:83-97.
%
% Inputs:
%   costMatrix - [nTracks x nDetections] cost matrix (Inf = forbidden)
%
% Outputs:
%   assignment - [nTracks x 1] detection index for each track (0 = unassigned)
%   cost       - Total assignment cost

[nRows, nCols] = size(costMatrix);

% Try MATLAB's built-in if available (Optimization Toolbox)
if exist('matchpairs', 'file') == 2
    % matchpairs available (R2019b+)
    maxCost = max(costMatrix(~isinf(costMatrix)));
    if isempty(maxCost), maxCost = 1; end
    C = costMatrix;
    C(isinf(C)) = maxCost * 100;  % Replace Inf with large cost
    
    try
        pairs = matchpairs(C, maxCost * 10);
        assignment = zeros(nRows, 1);
        cost = 0;
        for p = 1:size(pairs, 1)
            if ~isinf(costMatrix(pairs(p,1), pairs(p,2)))
                assignment(pairs(p,1)) = pairs(p,2);
                cost = cost + costMatrix(pairs(p,1), pairs(p,2));
            end
        end
        return;
    catch
        % Fall through to manual implementation
    end
end

% Manual implementation for environments without matchpairs
% Pad to square matrix
n = max(nRows, nCols);
C = zeros(n);
maxCost = max(costMatrix(~isinf(costMatrix)));
if isempty(maxCost), maxCost = 1; end
bigNum = maxCost * 1000;

C(1:nRows, 1:nCols) = costMatrix;
C(isinf(C)) = bigNum;
if nRows < n
    C(nRows+1:n, :) = bigNum;
end
if nCols < n
    C(:, nCols+1:n) = bigNum;
end

% Step 1: Subtract row minima
for i = 1:n
    C(i,:) = C(i,:) - min(C(i,:));
end

% Step 2: Subtract column minima
for j = 1:n
    C(:,j) = C(:,j) - min(C(:,j));
end

% Iterative assignment
maxIter = n * 10;
for iter = 1:maxIter
    % Find assignment via augmenting paths
    [rowAssign, colAssign] = find_assignment(C, n);
    
    nAssigned = sum(rowAssign > 0);
    if nAssigned >= n
        break;
    end
    
    % Find minimum uncovered element and adjust
    rowCovered = rowAssign > 0;
    colCovered = false(n, 1);
    
    % Mark columns covered by assigned rows
    for i = 1:n
        if rowAssign(i) > 0
            colCovered(rowAssign(i)) = true;
        end
    end
    
    % Find uncovered zeros and adjust
    uncoveredMin = inf;
    for i = 1:n
        for j = 1:n
            if ~rowCovered(i) && ~colCovered(j)
                uncoveredMin = min(uncoveredMin, C(i,j));
            end
        end
    end
    
    if isinf(uncoveredMin) || uncoveredMin == 0
        break;
    end
    
    % Subtract from uncovered, add to doubly covered
    for i = 1:n
        for j = 1:n
            if ~rowCovered(i) && ~colCovered(j)
                C(i,j) = C(i,j) - uncoveredMin;
            elseif rowCovered(i) && colCovered(j)
                C(i,j) = C(i,j) + uncoveredMin;
            end
        end
    end
end

% Extract assignment
assignment = zeros(nRows, 1);
cost = 0;
for i = 1:nRows
    j = rowAssign(i);
    if j > 0 && j <= nCols && ~isinf(costMatrix(i, j))
        assignment(i) = j;
        cost = cost + costMatrix(i, j);
    end
end

end


function [rowAssign, colAssign] = find_assignment(C, n)
% Greedy assignment on zero entries (simple heuristic)
rowAssign = zeros(n, 1);
colAssign = zeros(n, 1);
colUsed = false(n, 1);

% First pass: unique zeros in rows
for i = 1:n
    zeroIdx = find(C(i,:) < 1e-10 & ~colUsed');
    if numel(zeroIdx) == 1
        rowAssign(i) = zeroIdx;
        colAssign(zeroIdx) = i;
        colUsed(zeroIdx) = true;
    end
end

% Second pass: assign remaining
for i = 1:n
    if rowAssign(i) > 0, continue; end
    zeroIdx = find(C(i,:) < 1e-10 & ~colUsed');
    if ~isempty(zeroIdx)
        rowAssign(i) = zeroIdx(1);
        colAssign(zeroIdx(1)) = i;
        colUsed(zeroIdx(1)) = true;
    end
end

end
