function [optimalIndices, maxDiff_opt, bounds] = ...
    findIndividualOptimalTimePoints_interval_rSNR_slidingWindow( ...
    rSNR_left, rSNR_right, timePoints, subjectsForPlot, ...
    minlowerband, maxUpperband, varargin)
% FINDINDIVIDUALOPTIMALTIMEPOINTS_INTERVAL_RSNR_SLIDINGWINDOW
%
% For each subject, this function searches within the time range
% [minlowerband, maxUpperband] and finds a sub-interval of length
% 'windowSize' (in # of samples) that maximizes the area (sum) of
% |rSNR_left - rSNR_right|. 
%
% USAGE:
%   [optimalIndices, maxAUC, bounds] = findIndividualOptimalTimePoints_interval_rSNR_slidingWindow( ...
%       rSNR_left, rSNR_right, timePoints, [1,5], 0, 300, 'WindowSize', 10);
%
% INPUTS:
%   rSNR_left, rSNR_right : [N x T], N subjects, T time points
%   timePoints            : [T x 1] or [1 x T], time values for each column
%   subjectsForPlot       : vector of subject indices to plot (optional)
%   minlowerband, maxUpperband : numeric, bounding times (e.g., 0 and 300)
%
% OPTIONAL NAME-VALUE ARG:
%   'WindowSize' : (integer) number of consecutive samples for the window (default=5)
%   'UseTrapz'   : (logical) if true, uses trapz(...) for uneven time steps 
%                  instead of sum(...) (default=false).
%
% OUTPUTS:
%   optimalIndices : [N x 1], the *center* index (in original timePoints) of the best window
%   maxDiff_opt    : [N x 1], the maximum area (AUC) for each subject
%   bounds         : [N x 2], the [startIdx, endIdx] in original timePoints
%                    that define the best window for each subject
%
% -------------------------------------------------------------------------
% Author: ChatGPT (based on your original function)
% -------------------------------------------------------------------------

%% 1) Parse optional inputs
p = inputParser;
addParameter(p, 'WindowSize', 5, @(x) isnumeric(x) && x>0);
addParameter(p, 'UseTrapz',   false, @islogical);
parse(p, varargin{:});

windowSize = p.Results.WindowSize;
useTrapz   = p.Results.UseTrapz;

%% 2) Restrict to valid time range [minlowerband, maxUpperband]
% timePoints = timePoints(:)'; % force row vector
% validIdx   = find(timePoints >= minlowerband & timePoints <= maxUpperband);

validIdx = find(timePoints(:,1) >= minlowerband & timePoints(:,1) <= maxUpperband);

% Subset the rSNR arrays and the time vector
rSNR_left_valid  = rSNR_left(:,  validIdx);
rSNR_right_valid = rSNR_right(:, validIdx);
timeVec_valid    = timePoints(validIdx);

numSubjects  = size(rSNR_left, 1);
numValidTime = length(validIdx);

%% 3) Initialize outputs
optimalIndices = nan(numSubjects,1);  % center of best window in *original* index space
maxDiff_opt    = nan(numSubjects,1);  % the AUC value for that window
bounds         = nan(numSubjects,2);  % [startIndex, endIndex] in *original* space

%% 4) Loop over subjects
for subj = 1:numSubjects

    % Absolute difference in the valid range
    absDiff_valid = abs(rSNR_left_valid(subj,:) - rSNR_right_valid(subj,:));

    bestAUC    = -inf;
    bestStart  = NaN;
    bestEnd    = NaN;

    % Slide a window of size 'windowSize' over validIdx
    for startPos = 1:(numValidTime - windowSize + 1)
        endPos = startPos + windowSize - 1;

        % Extract the portion of the diff
        windowData = absDiff_valid(startPos:endPos);

        % Compute area
        if useTrapz
            % If time steps are NOT uniform, use numerical integration
            windowTime = timeVec_valid(startPos:endPos);
            areaVal    = trapz(windowTime, windowData);
        else
            % If time steps are uniform or you just want a sum
            areaVal    = sum(windowData);
        end

        if areaVal > bestAUC
            bestAUC   = areaVal;
            bestStart = startPos;
            bestEnd   = endPos;
        end
    end

    if ~isnan(bestStart)
        % best window found
        maxDiff_opt(subj)    = bestAUC;
        % map local indices back to the *original* timePoints indexing
        globalStart          = validIdx(bestStart);
        globalEnd            = validIdx(bestEnd);
        bounds(subj,:)       = [globalStart, globalEnd];

        % define "optimal index" as the approximate center of that window
        centerLocal          = floor((bestStart + bestEnd)/2);
        centerGlobal         = validIdx(centerLocal);
        optimalIndices(subj) = centerGlobal;
    end

end % end for subjects

%% 5) Optional plotting for selected subjects
if ~isempty(subjectsForPlot) && any(~isnan(subjectsForPlot))
    figure('Color','w','Position',[200,300,900,600]);
    for i = 1:length(subjectsForPlot)
        subj = subjectsForPlot(i);
        subplot(ceil(length(subjectsForPlot)/2), 2, i);

        % Full-series absolute difference
        absDiff_full = abs(rSNR_left(subj,:) - rSNR_right(subj,:));
        plot(timePoints, rSNR_left(subj,:),  'b-', 'DisplayName','rSNR Left'); hold on;
        plot(timePoints, rSNR_right(subj,:), 'r-', 'DisplayName','rSNR Right');
        plot(timePoints, absDiff_full,       'k-', 'LineWidth',1.5, 'DisplayName','|L-R|');

        % If we found a window
        if ~isnan(bounds(subj,1)) && ~isnan(bounds(subj,2))
            wStart = bounds(subj,1);
            wEnd   = bounds(subj,2);
            xx     = timePoints(wStart:wEnd);
            yy     = absDiff_full(wStart:wEnd);

            patch([xx, xx(end), xx(1)], [yy, 0, 0], [0.8 0.8 0.2], ...
                'FaceAlpha',0.3, 'EdgeColor','none', 'DisplayName','Optimal Window');

            title(sprintf('Subj %d | Window [%d:%d], AUC=%.1f', ...
                subj, wStart, wEnd, maxDiff_opt(subj)));
        else
            title(sprintf('Subj %d | No window found', subj));
        end

        xlabel('Time'); ylabel('rSNR / |L-R|');
        xlim([minlowerband, maxUpperband]);
        legend('Location','best');
        box off; set(gca,'Color','none');
    end
end

end
