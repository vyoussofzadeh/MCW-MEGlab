
function plotIndividualCounting(roi_cnt_pt, LI, ID, timePoints, subjectsForPlot)
% Plots active vertices count for selected subjects,
% and overlays the Laterality Index (LI) on the same plot.

% Ensure subjectsForPlot is within the ID array bounds
subjectsForPlot = subjectsForPlot(subjectsForPlot <= length(ID) & subjectsForPlot >= 1);

if isempty(subjectsForPlot)
    warning('No valid subjects specified for plotting.');
    return;
end

% Create a figure for all subplots
figure;

% Iterate over the specified subjects for plotting
for index = 1:length(subjectsForPlot)
    subj = subjectsForPlot(index);
    
    % Calculate the total number of subplots to determine the layout
    totalPlots = length(subjectsForPlot);
    numRows = ceil(totalPlots / 2);
    
    % Create subplot
    sp = subplot(numRows, 2, index); % Organize subplots in a 2-column layout
    
    yyaxis left; % Use left y-axis for roi_cnt_pt
    plot(timePoints, roi_cnt_pt(subj, :), 'LineWidth', 1.5); % Plot roi_cnt_pt as a solid line
    ylabel('Active Vertices Count');
    
    yyaxis right; % Use right y-axis for LI
    %     plot(timePoints, LI(subj, :), 'k', 'LineWidth', 1.5); % Overlay the LI as a dashed black line
    plot(timePoints, LI(subj, :), 'LineWidth', 1.5); % Overlay the LI as a dashed black line
    ylabel('Laterality Index (LI)');
    
    title(sprintf('Subject %s', ID{subj})); % Title with Subject ID
end
set(gcf, 'Position', [100   100   1200   1000]);

% Applying xlabel and legend to the last subplot used
% This is a workaround to apply these to the "figure" rather than a specific subplot
if exist('sp', 'var') % Check if 'sp' exists to avoid errors
    axes(sp); % Set the current axes to the last subplot for xlabel and legend application
    xlabel('Time Points');
    lgd = legend({'Active Vertices Count', 'LI'}, 'Location', 'bestoutside');
    lgdPos = lgd.Position; % Get current position
    lgdPos(2) = lgdPos(2) - 0.11; % Move legend down
    lgd.Position = lgdPos; % Set new position
end

end


