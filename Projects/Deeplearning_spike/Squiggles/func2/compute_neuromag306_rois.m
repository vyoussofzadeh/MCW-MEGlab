function [roi_idx, roi_names, lay, pos_xy] = compute_neuromag306_rois(MEG_data)
% compute_neuromag306_rois
%
% Derive simple lobar/hemispheric ROIs for Neuromag306 sensors based on the
% 2D layout (left/right × posterior?anterior stripes).
%
% INPUT:
%   MEG_data - FieldTrip data struct with .grad and .label fields
%
% OUTPUT:
%   roi_idx   - [nChan x 1] ROI index for each channel (1..Nroi)
%   roi_names - {1 x Nroi} cell array of ROI names:
%               {'L-Occ','R-Occ','L-Temp','R-Temp',...
%                'L-CentPar','R-CentPar','L-Front','R-Front'}
%   lay       - layout struct from ft_prepare_layout
%   pos_xy    - [nChan x 2] 2D positions used for the ROI assignment
%
% USAGE:
%   [roi_idx, roi_names] = compute_neuromag306_rois(MEG_data);
%   plot_meg_roi_heatmap(MEG_data, roi_idx, roi_names, true, 1);

    % --- 1. Get layout for Neuromag306 sensors in YOUR channel order ---
    cfg = [];
    cfg.grad    = MEG_data.grad;      % neuromag306
    cfg.channel = MEG_data.label;     % keep same order as MEG_data.trial
    lay = ft_prepare_layout(cfg);

    pos    = lay.pos;                 % [nChan x 2] 2D positions
    labels = lay.label;               % should match MEG_data.label (+ COMNT/SCALE)
    nChan  = numel(MEG_data.label);

    % keep only real MEG channels (skip COMNT/SCALE if present)
    labels = labels(1:nChan);
    pos    = pos(1:nChan, :);

    x = pos(:,1);
    y = pos(:,2);

    % --- 2. Define left/right and anterior/posterior stripes ---
    xmid = median(x);
    yq   = quantile(y, [0.25 0.5 0.75]);   % 4 stripes back?front

    % ROI order: 1=L-Occ, 2=R-Occ, 3=L-Temp, 4=R-Temp,
    %            5=L-CentPar, 6=R-CentPar, 7=L-Front, 8=R-Front
    roi_names = {...
        'L-Occ', 'R-Occ', ...
        'L-Temp','R-Temp', ...
        'L-CentPar','R-CentPar', ...
        'L-Front','R-Front'};

    roi_idx = zeros(nChan,1);  % ROI index per channel

    for ch = 1:nChan
        % stripe in y (posterior?anterior)
        if y(ch) < yq(1)
            stripe = 1;   % occipital
        elseif y(ch) < yq(2)
            stripe = 2;   % temporal
        elseif y(ch) < yq(3)
            stripe = 3;   % centralparietal
        else
            stripe = 4;   % frontal
        end

        % hemisphere in x
        if x(ch) < xmid
            hemi = 1;     % left
        else
            hemi = 2;     % right
        end

        roi_idx(ch) = (stripe-1)*2 + hemi;  % 1..8
    end

    if nargout > 3
        pos_xy = pos;
    end
end
