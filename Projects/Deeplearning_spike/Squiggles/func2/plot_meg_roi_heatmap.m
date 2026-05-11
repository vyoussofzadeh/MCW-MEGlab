function plot_meg_roi_heatmap(MEG_data, roi_idx, roi_names, do_zscore, trial_idx)
% plot_meg_roi_heatmap
%
% Plot a channel × time heatmap of MEG data, with channels grouped by ROI.
%
% INPUTS:
%   MEG_data  - FieldTrip data struct (continuous or epoched)
%   roi_idx   - [nChan x 1] ROI index for each MEG channel (1..Nroi)
%   roi_names - {1 x Nroi} cell array of ROI names
%   do_zscore - (optional) logical, z-score each channel over time (default: false)
%   trial_idx - (optional) which trial to plot (default: 1)
%
% EXAMPLE:
%   plot_meg_roi_heatmap(MEG_data, roi_idx, roi_names, true, 1);

    if nargin < 4 || isempty(do_zscore)
        do_zscore = false;
    end
    if nargin < 5 || isempty(trial_idx)
        trial_idx = 1;
    end

    % --- Extract data ---
    dat = MEG_data.trial{trial_idx};   % [nChan x nTime]
    t   = MEG_data.time{trial_idx};    % [1 x nTime]
    [nChan, ~] = size(dat);

    % --- Sanity checks ---
    if numel(roi_idx) ~= nChan
        error('roi_idx length (%d) does not match number of channels (%d).', ...
              numel(roi_idx), nChan);
    end

    if max(roi_idx) > numel(roi_names)
        error('roi_names has fewer entries (%d) than max(roi_idx) = %d.', ...
              numel(roi_names), max(roi_idx));
    end

    % --- Preprocess: z-score or abs ---
    if do_zscore
        dat_proc = zscore(dat, 0, 2);   % z-score along time
        ttl_tag  = ' (z-scored)';
    else
        dat_proc = dat;
        ttl_tag  = '';
    end

    % --- Reorder channels according to ROI index ---
    [roi_idx_sorted, sort_idx] = sort(roi_idx);
    dat_ord = dat_proc(sort_idx, :);

    % --- Plot heatmap ---
    figure;
    imagesc(t, 1:nChan, dat_ord);
    axis tight;
    set(gca, 'YDir', 'normal');       % so channel 1 is at bottom
    colorbar;
    xlabel('Time (s)');
    ylabel('Channels (grouped by ROI)');
    title(['MEG data channel × time, ROI-grouped', ttl_tag]);

    % --- Compute ROI boundaries and centers for labels ---
    edges   = [find(diff(roi_idx_sorted) ~= 0); numel(roi_idx_sorted)];
    starts  = [1; edges(1:end-1)+1];
    centers = floor((starts + edges)/2);

    yticks(centers);
    yticklabels(roi_names);
    set(gca, 'TickDir', 'out', 'FontSize', 10);

    hold on;
    for e = edges(1:end-1)'
        yline(e+0.5, 'k-');
    end
    hold off;
end
