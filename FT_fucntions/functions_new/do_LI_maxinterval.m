function Out = do_LI_maxinterval(cfg_main)


wi = cfg_main.wi;
LI_pt_new = cfg_main.LI_val;
net_sel = cfg_main.net_sel;
startTime = cfg_main.startTime;
endTime = cfg_main.endTime;

%%
timePoints = mean(wi,2);
% Define time range in seconds
startTime = startTime; % 200 ms
endTime = endTime; % 1 sec

% Find columns corresponding to the desired time range
startCol = find(timePoints >= startTime, 1, 'first');
endCol = find(timePoints <= endTime, 1, 'last');

% Determine the time point of max LI for each subject
LI_int = squeeze(LI_pt_new(net_sel,:,:));

%%
abs_LI = squeeze(abs(LI_int));
% abs_LI = squeeze((LI_int));

abs_LI2 = abs_LI(:,startCol:endCol);

[~, max_time_pts] = max(abs_LI2);

tsel = timePoints(startCol:endCol);

clc, close all
for i=1:size(abs_LI,1)
    
    %     figure, plot(tsel, abs_LI(i,startCol:endCol))
    %     xline(tsel(max_time_pts(i)),'r')
    figure, plot(abs_LI(i,:))
    hold on, xline((max_time_pts(i)) + startCol-1)
    
    figure, plot(timePoints, abs_LI(i,:))
    hold on, xline(tsel(max_time_pts(i)))
    pause,
end

%% Out
Out = [];
Out.max_time_pts = max_time_pts;
Out.tsel = tsel;
Out.max_time_idx = max_time_pts + startCol-1;

end

