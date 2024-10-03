function globalwindowthresh = do_calculateGlobalWindowThreshold(tmp, wi, idx_L, idx_R, cfg_main)
    % Calculate global window threshold based on amplitude data and HCP atlas regions
    %
    % Parameters:
    %   tmp       - Structure with fields 'ImageGridAmp' and 'Time'
    %   wi        - Window intervals as Nx2 matrix [start_time, end_time]
    %   idx_L     - Indices of left hemisphere regions
    %   idx_R     - Indices of right hemisphere regions
    %   cfg_main  - Configuration structure with field 'atlas'
    %
    % Returns:
    %   globalwindowthresh - Maximum threshold calculated across all specified windows

    % Look up indices for vertices or (sub)ROIs from HCP atlas
    if size(tmp.ImageGridAmp,1) > 360
        sScout = cfg_main.atlas;
        
        % Get left and right subregions from scout data
        LHscout = [];
        for i = 1:length(idx_L)
            LHscout = [LHscout, sScout.Scouts(idx_L(i)).Vertices];
        end
        
        RHscout = [];
        for i = 1:length(idx_R)
            RHscout = [RHscout, sScout.Scouts(idx_R(i)).Vertices];
        end
        idx_LR_updt = [LHscout, RHscout];
    else
        idx_LR_updt = [idx_L, idx_R];
    end

    % Calculate mean across windows
    mdwin = [];
    for j = 1:size(wi,1)
        timind1 = nearest(tmp.Time, wi(j,1));
        timind2 = nearest(tmp.Time, wi(j,2));
        dwin = tmp.ImageGridAmp(:,timind1:timind2);
        mdwin(j,:) = mean(dwin(idx_LR_updt,:),2);
    end

    % Calculate global threshold as the maximum mean amplitude
    globalwindowthresh = max(mdwin(:));

    % Optionally plot window mean values
    % figure, plot(wi(:,1),mdwin);
end

% Helper function to find nearest index of a time value
function idx = nearest(timeArray, timeValue)
    [~, idx] = min(abs(timeArray - timeValue));
end
