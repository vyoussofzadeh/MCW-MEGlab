
savepath = fullfile(outd.sub,'speech',['speech_',subj,'.mat']);
if exist(savepath, 'file') == 2
    load(savepath)
    disp('speech data was loaded')
else
    
    %% Reading event files
    cfg                         = [];
    cfg.dataset                 = datafile;
    raw_data = ft_preprocessing(cfg);
    label = raw_data.label;
    idx = find(strcmp(label,'MISC001'));
    if ~isempty(idx)
        
        %-
        load(['f_',subj,'.mat']);
        
        cfg                         = [];
        cfg.dataset                 = datafile;
        cfg.trialfun                = 'ft_trialfun_general'; % this is the default
        cfg.trialdef.eventtype      = f_data.cfg.trialdef.eventtype;
        cfg.trialdef.eventvalue     = f_data.cfg.trialdef.eventvalue; % the value of the stimulus trigger for fully incongruent (FIC).
        cfg.trialdef.prestim        = 1; % in seconds
        cfg.trialdef.poststim       = 3; % in seconds
        cfg = ft_definetrial(cfg);
        
        cfg.channel = {'MISC001'};
        cfg.hpfreq = 70;
        cfg.demean = 'yes';
        speech_data = ft_preprocessing(cfg);
        
        
        %%
        % %
        %     %
        %     cfg = [];
        %     cfg.viewmode = 'vertical';
        % %     cfg.continuous = 'no';
        % %     cfg.trials     = report.btrl;
        %     cfg.channel   = {'STI*'};
        %     ft_databrowser(cfg,raw_data);
        
        %%
        % sampleinfo =
        % cfg = [];
        % cfg.resamplefs = 1000;
        % speech_data = ft_resampledata(cfg, speech_data);
        
        textfile_rej = 'speech/bad_speech';
        %- Speech onset detection
        tt=[];
        for i=1:length(speech_data.trial)
            
            tmp = speech_data.trial{1,i} - mean(speech_data.trial{1,i});
            tmp = detrend(tmp);
            tmp1 = zeros(1,length(tmp));
            tmp1(1,200:end-200) = tmp(1,200:end-200);
            
            %     [thres_buf,env, bin] = envelop_hilbert_modified(abs(tmp1),20,1,20,1);
            [thres_buf,env, bin] = envelop_hilbert_modified(abs(tmp1),20,1,20,0);
            
            [a,b] = find(thres_buf > 0.8.*max(thres_buf));
            [d,initCross,finalCross,nextCross,midRef] =  dutycycle(bin());
            if isempty(initCross)
                idx = find(bin > 0); ipoints = idx(1);
            else
                max_idx = intersect(find(initCross < b(1)),find(finalCross > b(1)));
                if isempty(max_idx)
                    max_idx = find(initCross < b(1));
                    max_idx = max_idx(end);
                end
                ipoints = round(initCross(max_idx));
            end
            tt(i) = speech_data.time{1}(ipoints);
            ipoints_all(i) = ipoints;
            ipoints_good(i) = ipoints > 0.5.*mean(ipoints_all);
            
            figure, plot(speech_data.time{1}, abs(hilbert(speech_data.trial{1,i}))),
            hold on,
            scaled = (ipoints * (abs(min(speech_data.time{1,i})) + abs(max(speech_data.time{1}))))/length(speech_data.trial{1,i});
            scaled = scaled - abs(min(speech_data.time{1,i}));
            vline(scaled,'g',['speech onset:', num2str(scaled),'sec']),
            box off;
            title(['Trl ',num2str(i)])
            set(gca,'color','none');
            
            %     disp('1: OK');
            %     disp('2: A better guess');
            %     disp('3: Bad')
            %     onset_ok = input('how is the onset?');
            %     onset_ok = 2;
            %     switch onset_ok
            %         case 1
            %             good_onset(i)=1;
            %         case 2
            disp('Enter B if trial is bad, and G for good')
            [ot,~, onset_ok] = ginput(1);
            vline(ot,'g',['Onset guess:', num2str(scaled),'sec']),
            switch onset_ok
                case 103 % good
                    good_onset(i)=1;
                    tt(i) = tt(i);
                case 98 % bad
                    good_onset(i)=0;
                    tt(i) = tt(i);
                otherwise %man selection
                    good_onset(i)=1;
                    tt(i) = ot;
            end
            
            %             pause,
            %             tt(i) = ot;
            %         case 3
            %             good_onset(i) = 0;
            %     end
            close all,
        end
        
        idx_good_2 = find(good_onset==1);
        % tt1 = tt(idx);
        % ipoints_all1 = ipoints_all(idx);
        
        %%
        thre  = speech_data.time{1}(round(0.5.*mean(ipoints_all)));
        if thre < 0.4, thre = 0.4; end
        
        L = length(tt);
        idx_good = find(tt>= 0.5.*mean(tt));
        idx_good_1 = find(tt >= 0.4); % responses after 400ms are considered valid.
        idx_good = intersect(idx_good,idx_good_1);idx_good = intersect(idx_good,idx_good_2);
        idx_bad  = find(tt < 0.5.*mean(tt));
        
        figure,
        plot(tt, '*'),
        hold on
        tt1 = tt; tt1(idx_good) = nan;
        plot(tt1, 'r*'),
        hline(thre,'b',['mean:', num2str(thre),'sec']),
        box off;
        set(gca,'color','none');
        title('Speech onset'),
        ylabel('Time (sec)');
        xlabel('Trials');
        set(gca,'Xtick', 1:L,'XtickLabel',1:L);
        xlim([1 L]);
        set(gca,'FontSize',10,'XTickLabelRotation',90);
        grid
        
        %%
        samplepoints = tt.*speech_data.fsample;
        
        %%
        print('speech/speech','-dpng');
        badSpeech = table(idx_bad');
        badSpeech.Properties.VariableNames{'Var1'} = 'bSpeechs';
        writetable(badSpeech,textfile_rej,'Delimiter',' ');
        
        save(savepath,'tt','idx_good_2','speech_data')
        
        
    else
        disp('=====');
        disp('no speech data are available');
    end
end
%%
[C,ia,ib] = intersect(cln_data.sampleinfo(:,1),speech_data.sampleinfo(idx_good,1));

cfg = [];
cfg.trials = ia;
cln_data = ft_preprocessing(cfg, cln_data);

k=1;
nfig = floor(length(speech_data.trial)/12);
for i=1:nfig
    figure,
    for j=1:12
        subplot(3,4,j)
        plot(speech_data.time{1}, abs(hilbert(speech_data.trial{1,k}))),
        hold on,
        vline(tt(k),'r',['speech onset:', num2str(tt(k)),'sec']),
        if ~~find(idx_good_2==k)
            title(num2str(k))
        else
            title('BAD')
        end
        k=k+1;
    end
end



%%
% print('ica/ica2','-dpng');


% tt=[];
% for i=1:length(speech_data.trial)
%
%     tmp = speech_data.trial{1,i} - mean(speech_data.trial{1,i});
%
%     tmp = detrend(tmp);
% %     tmp = tmp - detrend_sdata;
%
%     tmp1 = zeros(1,length(tmp));
%     tmp1(1,200:end-200) = tmp(1,200:end-200);
%
% %     [mx, idx] = max(tmp);
%     [thres_buf,env, bin] = envelop_hilbert_modified(abs(tmp1));
%
%     [a,b] = find(thres_buf > 0.8.*max(thres_buf));
% %     hold on
% % %     plot(b,a,'*')
% %     y = ylim; % current y-axis limits
% %     plot([b(1) b(1)],[y(1) y(2)])
%
%     [d,initCross,finalCross,nextCross,midRef] =  dutycycle(bin);
%     if isempty(initCross)
%         idx = find(bin > 0); ipoints = idx(1);
%     else
%         max_idx = intersect(find(initCross < b(1)),find(finalCross > b(1)));
%         if isempty(max_idx)
%             max_idx = find(initCross < b(1));
%             max_idx = max_idx(end);
%         end
%         ipoints = round(initCross(max_idx));
%     end
%
% %
% %     [mx,idx] = max(d);
%
%     %     [ipoints, residual] = findchangepts(abs(hilbert(thres_buf)));
%     %     [ipoints, residual] = findchangepts((thres_buf));
% %     [ipoints, residual] = findchangepts(thres_buf);
%
%     tt(i) = speech_data.time{1}(ipoints);
%     ipoints_all(i) = ipoints;
%         ipoints_good(i) = ipoints > 0.5.*mean(ipoints_all);
%
%
%     figure, plot(speech_data.time{1}, abs(hilbert(speech_data.trial{1,i}))),
%     %             figure, plot(speech_data.time{1}, envelope(speech_data.trial{1,i})),
%     hold on,
%     scaled = (ipoints * (abs(min(speech_data.time{1,i})) + abs(max(speech_data.time{1}))))/length(speech_data.trial{1,i});
%     scaled = scaled - abs(min(speech_data.time{1,i}));
%     vline(scaled,'g',['speech onset:', num2str(scaled),'sec']),
%     box off;
%     set(gca,'color','none');
% %                 pause
%
% end