function plot_roi_rms_time(MEG_data, roi_idx, roi_names, trial_idx, do_zscore)
% plot_roi_rms_time
%
% Compute and plot RMS of MEG data within each ROI over time.
%
% INPUTS:
%   MEG_data  - FieldTrip data struct
%   roi_idx   - [nChan x 1] ROI index for each channel (1..Nroi)
%   roi_names - {1 x Nroi} cell array of ROI names
%   trial_idx - (optional) trial index to use (default: 1)
%   do_zscore - (optional) z-score each ROI time series across time (default: false)
%
% OUTPUT:
%   (none)  makes a figure: ROI × time heatmap of RMS amplitude
%
% EXAMPLE:
%   plot_roi_rms_time(MEG_data, roi_idx, roi_names);
%
%   plot_roi_rms_time(MEG_data, roi_idx, roi_names, 1, true);  % z-scored

    if nargin < 4 || isempty(trial_idx)
        trial_idx = 1;
    end
    if nargin < 5 || isempty(do_zscore)
        do_zscore = false;
    end

    dat = MEG_data.trial{trial_idx};     % [nChan x nTime]
    t   = MEG_data.time{trial_idx};      % [1 x nTime]
    [nChan, nTime] = size(dat);

    % --- sanity checks ---
    if numel(roi_idx) ~= nChan
        error('roi_idx length (%d) does not match number of channels (%d).', ...
              numel(roi_idx), nChan);
    end

    Nroi = numel(roi_names);
    if max(roi_idx) > Nroi
        error('roi_names has fewer entries (%d) than max(roi_idx) = %d.', ...
              Nroi, max(roi_idx));
    end

    % --- compute RMS per ROI over channels ---
    roi_rms = nan(Nroi, nTime);   % [Nroi x nTime]

    for r = 1:Nroi
        ch_idx = find(roi_idx == r);
        if isempty(ch_idx)
            continue;
        end
        % RMS over channels within this ROI at each time point
        roi_dat = dat(ch_idx, :);             % [nCh_roi x nTime]
        roi_rms(r, :) = sqrt(mean(roi_dat.^2, 1));
    end

    % --- optional z-score per ROI across time ---
    if do_zscore
        roi_rms = zscore(roi_rms, 0, 2);
        ttl_tag = ' (z-scored)';
    else
        ttl_tag = '';
    end

    % --- plot ROI x time heatmap ---
    figure;
    imagesc(t, 1:Nroi, roi_rms);
    axis tight;
    colorbar;
    xlabel('Time (s)');
    ylabel('ROI');
    set(gca, 'YDir', 'normal');
    set(gca, 'TickDir', 'out', 'FontSize', 10);
    yticks(1:Nroi);
    yticklabels(roi_names);
    title(['ROI RMS over time', ttl_tag]);
end
