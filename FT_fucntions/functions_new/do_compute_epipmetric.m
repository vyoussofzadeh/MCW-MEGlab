function [time_occur, val] = do_compute_epipmetric(cfg, data)


LL = length(data);

%% 
if ~isempty(cfg.foi)
    [~, mn1] = min(abs(data{1}.freq - cfg.foi(1)));
    [~, mn2] = min(abs(data{1}.freq - cfg.foi(2)));
    for i = 1:LL
        data{i}.pow = data{i}.pow(:,mn1:mn2);
        data{i}.freq = data{i}.freq(mn1:mn2);
    end
end

%%
t_info_all = [];
figure,
for i = 1:LL
    subplot(LL,1,i)
    
    switch cfg.metric
        case 'pow'
            val = mean(data{i}.pow,2);
            val = val./max(val);
            plot(data{1}.time,val);
        case 'kurt'
            val = kurtosis(data{i}.pow, [], 2);
            plot(data{1}.time,val);
        case 'skewness'
            val = skewness(data{i}.pow, [], 2);
            plot(data{1}.time,val);
    end
    title(cfg.labels{i})
    
    %%
%     TF = islocalmax(val,'MinProminence',2);
    TF = islocalmax(val);
    [mx, ~] = max(TF);
    idd  = find(TF  == mx);
    t_info = data{1}.time(idd);
    

%     [PKS, peaksample] = findpeaks(val, 'MinPeakDistance', 500); % peaks have to be separated by 300 sample points to be treated as separate
%     
%     [mx, ~] = max(val);
%     idd  = find(val  > 0.8.*mx);
%     t_info = data{1}.time(idd);
    
    %%
    %     disp([cfg.labels{i},':']);
    %     disp(t_info)
    t_info_all = [t_info_all, t_info];
    
    
    %%
    for j=1:length(t_info)
        hold on
        plot([t_info(j) t_info(j)],[min(val) max(val)],'-', 'Color',[0.5 0.5 0.5])
    end
end
if LL > 5, set(gcf, 'Position', [1110   300   300   1200]); end

if length(t_info_all) > 2
    t_info_all = sort(t_info_all, 'ascend');
    diff(t_info_all);
    idx = diff(t_info_all) > 3;
    time_occur  = t_info_all(idx);
else
    time_occur  = t_info_all;
end
disp('identified peak:')
disp(time_occur);


end