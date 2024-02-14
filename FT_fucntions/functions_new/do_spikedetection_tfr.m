function [time_occur, out_val] = do_spikedetection_tfr(cfg, val)

fs = cfg.fsample;

%%
% sens = sensorlist();

%%
% disp([cfg.art.eog./fs])
% disp([cfg.art.jump./fs])
for i=1:size(cfg.art.eog,1)
    [~, idx1] = min(abs(val.time - cfg.art.eog(i,1)/fs));
    [~, idx2] = min(abs(val.time - cfg.art.eog(i,2)/fs));
    val.pow(:,idx1:idx2) = nan;
end
for i=1:size(cfg.art.jump,1)
    [~, idx1] = min(abs(val.time - cfg.art.jump(i,1)/fs));
    [~, idx2] = min(abs(val.time - cfg.art.jump(i,2)/fs));
    val.pow(:,idx1:idx2) = nan;
end
for i=1:size(cfg.art.rejseg,1)
    [~, idx1] = min(abs(val.time - cfg.art.rejseg(i,1)/fs));
    [~, idx2] = min(abs(val.time - cfg.art.rejseg(i,2)/fs));
    val.pow(:,idx1:idx2) = nan;
end

%%
[~, minidx1] = min(abs(val.freq - cfg.foi(1))); val.freq(minidx1);
[~, minidx2] = min(abs(val.freq - cfg.foi(2))); val.freq(minidx2);


ts_pp = (abs(hilbert(abs(val.pow(minidx1:minidx2,:)))));
ts_pp = ts_pp./max(ts_pp(:));


d_in = nanmean(ts_pp,1);
% d_in = nanmedian(ts_pp,1);
% d_in = nanstd(ts_pp,1);

% d_in = d_in - mean(d_in);

%%
thre_ok = 0;
while thre_ok == 0
    disp('enter the threshold level (e.g, 0.2):')
    inp.thre = input('');
    thre = inp.thre;
    [mx, ~] = max(d_in(3:end));
    peaksample  = d_in  > thre.*mx;
    time_occur = val.time(peaksample);
    out_val = d_in(peaksample);
    
    % k = 1;
    % time_occur_sel = [];
    % out_val_sel = [];
    % if length(time_occur) > 2
    %     time_occur = sort(time_occur, 'ascend');
    %     for i= 1:length(time_occur)-1
    %         df = time_occur(i+1) - time_occur(i);
    %         if df > 5
    %             time_occur_sel(k) = time_occur(i);
    %             out_val_sel(k) = out_val(i);
    %             k = k + 1;
    %         end
    %     end
    % else
    %     time_occur_sel  = time_occur;
    % end
    %
    % disp(time_occur_sel)
    % time_occur = time_occur_sel;
    % out_val = out_val_sel;
    
    %         diff(time_occur);
    %         idx = find(diff(time_occur) > 5);
    %         time_occur  = time_occur(idx+1);
    %         out_val = out_val(idx+1);
    
    %%
    % [~, peaksample] = findpeaks(d_in, 'MinPeakDistance', 10); % peaks have to be separated by 300 sample points to be treated as separate
    % time_occur = val.time(peaksample);
    % out_val = d_in(peaksample);
    % idx = out_val > 0.3.*max(out_val);
    % time_occur = time_occur(idx);
    % out_val = out_val(idx);
    
    %%
    if cfg.plot == 1
        
        figure,plot(val.time, d_in)
        for j=1:length(time_occur)
            hold on
            plot([time_occur(j) time_occur(j)],[min(d_in), max(d_in)],'-', 'Color',[0.7 0.7 0.7]);
            text(time_occur(j),out_val(j), num2str(j))
        end
        plot([val.time(1) val.time(end)],[thre.*mx, thre.*mx],'-', 'Color',[0.8 0.8 0.1]);
        title([cfg.ttl, ' (amp envelope), thre: ', num2str(thre),'.*max'])
        set(gcf, 'Position', [800   500   1000   200]);
        
        %     figure,
        %     subplot 211
        %     plot(val.time, mean(ts_pp,1)), xlabel('Time'),ylabel('mean pow');
        %     subplot 212
        %     plot(val.freq(minidx1:minidx2), mean(ts_pp,2)), xlabel('Freq'), ylabel('mean pow');
        
        kk =  (cfg.foi(2) - cfg.foi(1))/10;
        
        figure, imagesc(val.time, val.freq(minidx1:minidx2), ts_pp);
        ylabel('Freq (Hz)'); xlabel('Time (sec)')
        axis xy;
        %     for j=1:length(time_occur)
        %         hold on
        %         plot([time_occur(j) time_occur(j)],[cfg.foi(2)-kk, cfg.foi(2)],'-', 'Color',[1 1 1]);
        %         plot([time_occur(j) time_occur(j)],[cfg.foi(1), cfg.foi(1)+kk],'-', 'Color',[1 1 1]);
        %         text(time_occur(j),out_val(j), num2str(j))
        %     end
        title([cfg.ttl, ' (TFR), thre: ', num2str(thre),'.*max'])
        set(gcf, 'Position', [800   800   1000   200]);
    end
    % disp(time_occur)
    % disp('====')
    
    
    disp('OK with the threshold (yes:y, no:n)?')
    ask_thre_ok = input('','s');
    switch ask_thre_ok
        case 'y'
            thre_ok = 1;
        otherwise
            close all
    end
end
end
