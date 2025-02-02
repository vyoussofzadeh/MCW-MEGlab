function plotMEGLITimeseriesWithRSNR(MEG_LI_matrix, rSNR_left, rSNR_right, timeVec, discordSubs, subIDs)
% PLOTMEGLITIMESERIESWITHRSNR Plots the per-subject time series of MEG LI and
% (rSNR_left - rSNR_right) for a list of discordant subjects.
%
%   plotMEGLITimeseriesWithRSNR(MEG_LI_matrix, rSNR_left, rSNR_right, timeVec, discordSubs, subIDs)
%
%   INPUTS:
%       MEG_LI_matrix - (#Subjects x #TimePoints) numeric array of MEG LI values
%                       over time for each subject.
%       rSNR_left     - (#Subjects x #TimePoints) numeric array (left hemisphere rSNR).
%       rSNR_right    - (#Subjects x #TimePoints) numeric array (right hemisphere rSNR).
%       timeVec       - (#TimePoints x 1) numeric vector of time points (e.g. midpoints),
%                       must match columns in the above arrays.
%       discordSubs   - Vector of subject indices that are considered "discordant".
%       subIDs        - (#Subjects x 1) numeric/cell array of subject labels (or 1:N).
%
%   This function:
%       1) Creates a figure.
%       2) For each subject in discordSubs, generates a subplot:
%           (a) Plots the MEG LI over time in one color.
%           (b) Overlays (rSNR_left - rSNR_right) in another color.
%           (c) Adds subject label, axis labels, and a small legend.
%
%   EXAMPLE:
%       plotMEGLITimeseriesWithRSNR(MEG_LI_dyn, rSNR_L, rSNR_R, wi_mid, [4,5,10], 1:72);
%
%   Author: (Your Name / Organization)
%   Date: (Date)

    if nargin < 6
        error('All 6 inputs are required. Check your usage.');
    end
    
    % Check dimensions
    [numSubjects, numTimePoints] = size(MEG_LI_matrix);
    if size(rSNR_left,1) ~= numSubjects || size(rSNR_left,2) ~= numTimePoints
        error('rSNR_left must match (#Subjects x #TimePoints) of MEG_LI_matrix.');
    end
    if size(rSNR_right,1) ~= numSubjects || size(rSNR_right,2) ~= numTimePoints
        error('rSNR_right must match (#Subjects x #TimePoints) of MEG_LI_matrix.');
    end
    if length(timeVec) ~= numTimePoints
        error('timeVec length must match the #TimePoints of MEG_LI_matrix.');
    end
    if length(subIDs) ~= numSubjects
        error('subIDs must have the same length as #Subjects.');
    end

    % Create a new figure
    figure('Name','MEG LI + rSNR difference (Discordant Subjects)','Color','w');
    nDiscord = length(discordSubs);

    % Decide subplot arrangement
    % e.g., row x col = (ceil(nDiscord/2) x 2)
    nRows = ceil(nDiscord/2);
    nCols = 2;

    for i = 1:nDiscord
        subIdx = discordSubs(i);  % subject index in 1..numSubjects

        % Subplot
        subplot(nRows, nCols, i);
        hold on; box off; grid on;

        % Plot MEG LI over time (blue line)
        plot(timeVec, MEG_LI_matrix(subIdx,:), 'b-', 'LineWidth',1.5, ...
            'DisplayName','MEG LI');

        % Plot rSNR difference (left - right) in red
        rsnrDiff = rSNR_left(subIdx,:) - rSNR_right(subIdx,:);
        plot(timeVec, rsnrDiff, 'r--', 'LineWidth',1.5, ...
            'DisplayName','rSNR (L-R)');

        % Add subject label in the title
        thisSubLabel = num2str(subIDs(subIdx));  % convert to string if numeric
        title(sprintf('Subject %s (Discordant)', thisSubLabel), 'FontSize',10);

        % Axis labels (optional)
        xlabel('Time');
        ylabel('Value');

        % Optional reference lines
        yline(0,'k-','HandleVisibility','off');
        % xline(...) if you have an event or time=0 reference

        legend('Location','best','Box','off');
    end

    % Overall styling
    set(gcf, 'Position',[200,200,1200,600]);
    % axis tight per subplot -> automatically done by MATLAB
end
