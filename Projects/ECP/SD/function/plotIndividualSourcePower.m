function plotIndividualSourcePower(powerLeft, powerRight, MEG_LI, fMRI_LI, ID, timePoints, LI_opt,subjectsForPlot)
% Plots power values for left and right hemispheres for selected subjects,
% and overlays the Laterality Index (LI) on the same plot.

% Ensure subjectsForPlot is within bounds
subjectsForPlot = subjectsForPlot(subjectsForPlot <= length(ID) & subjectsForPlot >= 1);

numSubjects = length(subjectsForPlot);

% Create a figure outside the loop for all subplots if there are subjects to plot
if ~isempty(subjectsForPlot) && any(~isnan(subjectsForPlot))
    figure;
end

% Iterate over the specified subjects for plotting
for index = 1:numSubjects
    
    subj = subjectsForPlot(index);
    
    % Optional Plotting for this subject
    subplot(ceil(numSubjects / 2), 2, index); % Organize subplots in a 2-column layout
    yyaxis left; % Use left y-axis for power
    plot(timePoints, powerLeft(subj, :), 'b'); % Plot left hemisphere power as a blue line
    hold on;
    plot(timePoints, powerRight(subj, :), '-r'); % Plot right hemisphere power as a red line
    ylabel({'dB'});
    
    yyaxis right; % Use right y-axis for LI
    plot(timePoints, MEG_LI(subj, :), 'k', 'LineWidth', 1.5); % Overlay the LI as a dashed black line
    ylabel('LI');
    
    % Automatically adjust the x-axis based on the full range of timePoints
    % The line setting xlim([lowerBound upperBound]); has been removed
    
    hold off;
    title(sprintf('Subject %s', ID{subj}));
    title(sprintf('Subject %s', num2str(subj)));
    title(sprintf('S %d | fMRI LI: %.2f | MEG LI: %.2f', subj, fMRI_LI(subj), LI_opt(index)));
    
end

% Move xlabel and legend outside the loop to apply to the last subplot
% Set common xlabel for the figure
xlabel('Time Points');
% Place the legend in a consistent location
lgd = legend({'Left Power', 'Right Power', 'LI'}, 'Location', 'southoutside');
lgdPos = lgd.Position; % Get current position
lgdPos(2) = lgdPos(2) - 0.11; % Move legend down
lgd.Position = lgdPos; % Set new position

end

