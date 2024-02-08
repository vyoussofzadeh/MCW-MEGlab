function Out = do_LI_maxinterval(cfg_main)

% Extracting configuration settings
wi = cfg_main.wi;
LI_pt_new = cfg_main.LI_val;
net_sel = cfg_main.net_sel;
startTime = cfg_main.startTime;
endTime = cfg_main.endTime;

% Average time points across intervals
timePoints = mean(wi, 2);

% Extract Lateralization Index for selected network
LI_int = squeeze(LI_pt_new(net_sel, :, :));

% Define the selected interval
startCol = find(timePoints >= startTime, 1, 'first');
endCol = find(timePoints <= endTime, 1, 'last');
LI_int_sel = LI_int(:, startCol:endCol);

% Find the time point of maximum LI within the selected interval
[~, max_time_pts_sel] = max(abs(LI_int_sel), [], 2);
% [~, max_time_pts_sel] = max((LI_int_sel), [], 2);
max_time_vals_sel = timePoints(startCol - 1 + max_time_pts_sel);

% Adjusting max_time_pts_sel to reflect their positions in the full LI array
max_time_pts_full_range = max_time_pts_sel + startCol - 1;

% Output
Out = struct();
Out.max_time_pts_sel = max_time_pts_full_range; % Indices of max time points in full time range
Out.max_time_vals_sel = max_time_vals_sel;      % Time stamps of max LI points in selected interval

% Sanity Check Plot
% for i = 1:size(LI_int, 1)
%     figure; hold on;
%     plot(timePoints, LI_int(i, :));
%     xline(Out.max_time_vals_sel(i), 'r', 'LineWidth', 1.5);
%     xlabel('Time (s)');
%     ylabel('Lateralization Index');
%     title(['Sanity Check: Max LI Point in Selected Interval for Subject ' num2str(i)]);
%     legend('LI Time Course', 'Max LI Point', 'Location', 'best');
%     hold off;
%     pause;  % Pause to view each figure
% end

end
