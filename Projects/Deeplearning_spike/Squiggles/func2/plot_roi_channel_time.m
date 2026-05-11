function plot_roi_channel_time(data_mat, time_vec, roi_idx, roi_names, fig_title, dt_tick)
% plot_roi_channel_time
%
% Plot a channels × time matrix, grouped by ROIs on the y-axis.
%
% INPUTS:
%   data_mat  - [nChan x nTime] matrix (e.g. MEG, spike prob, etc.)
%   time_vec  - [1 x nTime] time vector (s)
%   roi_idx   - [nChan x 1] ROI index for each channel (1..Nroi)
%   roi_names - {1 x Nroi} cell array of ROI names
%   fig_title - (optional) figure title (default:
%               'Channel x Time grouped by Neuromag306 ROIs')

if nargin < 5 || isempty(fig_title)
    fig_title = 'Channel x Time grouped by Neuromag306 ROIs';
end

[nChan, nTime] = size(data_mat);

% --- sanity checks ---
if numel(time_vec) ~= nTime
    error('Length of time_vec (%d) does not match nTime (%d).', ...
        numel(time_vec), nTime);
end
if numel(roi_idx) ~= nChan
    error('roi_idx length (%d) does not match number of channels (%d).', ...
        numel(roi_idx), nChan);
end
if max(roi_idx) > numel(roi_names)
    error('roi_names has fewer entries (%d) than max(roi_idx) = %d.', ...
        numel(roi_names), max(roi_idx));
end

% --- sort by ROI and reorder data ---
[roi_idx_sorted, sort_idx] = sort(roi_idx);
out_ord = data_mat(sort_idx, :);

% --- plot ---
% figure;

% --- plot ---
figure('Units','normalized','Position',[0.2 0.4 0.6 0.3]); 
% [left bottom width height] in normalized units (01)

imagesc(time_vec, 1:nChan, out_ord);
axis tight;
% colorbar;
xlabel('Time (s)');
ylabel('Channels (grouped by ROI)');
title(fig_title);
set(gca, 'YDir', 'normal');   % channel 1 at bottom
set(gca, 'TickDir', 'out', 'FontSize', 9);

% --- ROI boundaries & labels ---
edges   = [find(diff(roi_idx_sorted) ~= 0); numel(roi_idx_sorted)];
starts  = [1; edges(1:end-1) + 1];
centers = floor((starts + edges) / 2);

yticks(centers);
yticklabels(roi_names);

hold on;
for e = edges(1:end-1)'
    yline(e + 0.5, 'w-');
end
hold off;

colormap(gca, 'jet');      % or just colormap('hot')
% colorbar;
colorbar off

%% >>> finer x-axis legend values <<<
Tstart = time_vec(1);
Tend   = time_vec(end);

% choose spacing (e.g., every 1 s; change to 0.5 for 500 ms, etc.)
% dt_tick = 5;   % seconds

if ~isempty(dt_tick)
    xticks(Tstart:dt_tick:Tend);
    xtickformat('%.1f');    % show 1 decimal (e.g., 0.0, 1.0, 2.0, ...)
end
%% <<< end x-axis control <<<


end
