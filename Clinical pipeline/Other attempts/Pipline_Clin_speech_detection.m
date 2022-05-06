clear; clc, close('all'); warning off

%% Flags
flag.preprocessing.filtering = 1;
flag.preprocessing.artifact = 1;
flag.preprocessing.ica = 1;
flag.notch = 1;
flag.freq = 1;     % TFR & FFT
flag.time = 1;     % Time-locked & Cov estimation
flag.gave = 0;     % grand average analysis
flag.anatomy = 1;     % grand average analysis
flag.sourceanalysis = 1;     % grand average analysis
flag.speechanalysis = 1;     % speech analysis

%% Initial settings
set(0,'DefaultFigureWindowStyle','docked')

cd '/MEG_data/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));

%- Input dir
indir = '/MEG_data/epilepsy';
%- Output dir
outdir = '/MEG_data/Vahab/Processed_data';

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/MEG_data/Vahab/Github/tools';
[allpath, atlas] = vy_init(cfg_init);

%%
disp('1: Definition naming')
disp('2: Picture naming');
disp('3: Motor');
task = input('Eneter the task: ');
switch task
    case 1
        %- Auditory definition naming
        tag = 'DFN';
    case 2
        %- Visual picture naming
        tag = 'PN';
    case 3
        %- Motor task
        tag = 'motor';
end

%%
cd(indir)
[subjdir] = uigetdir;

%%
d = rdir([subjdir,['/**/','sss','/*',tag,'*/*raw_tsss.fif']]);

%%
clear subj datafolder datafile datafile1
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    datafolder{i} = pathstr;
    datafile{i} = d(i).name;
    Index = strfind(datafile{i}, '/');
    subj = datafile{i}(Index(3)+1:Index(4)-1);
end
datafile1 = datafile';
disp(datafile1)
if length(datafile1) > 1 
   datasel = input('choose which data to analyze, row number:'); 
else
    datasel = 1;
end
disp([subj, ' and,'])
disp([datafile1{datasel}, 'was selected for the analysis ...'])
disp('============');

%%
epoch_type = 'STI101';

%% 4D layout
cfg = [];
cfg.layout = 'neuromag306mag.lay';
lay = ft_prepare_layout(cfg);
% ft_layoutplot(cfg);
disp('============');

%%
close all
datafile = datafile1{datasel}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
Index = strfind(datafile, '/');
Date  = datafile(Index(4)+1:Index(5)-1);
disp('============');
disp(datafile)
disp(['subj:',subj])
disp(['Date:',Date])
disp('============');

%%
%-elec/grad
sens = ft_read_sens(datafile);
sens = ft_convert_units(sens,'mm');

%%
outd.sub = fullfile(outdir,'ft_process',subj, tag);
if exist(outd.sub, 'file') == 0
    mkdir(outd.sub);   %create the directory
end
cd(outd.sub)
disp(['outputdir:',outd.sub])
disp('============');

%% Speech
savepath = ('speech');
if exist(savepath, 'file') == 0, mkdir(savepath), end

%%
% load(['f_',subj,'.mat']);

hdr = ft_read_header(datafile);
Index = strfind(hdr.label,{'MISC001'});
Index = find(not(cellfun('isempty',Index)));

%%
event = ft_read_event(datafile);
clear val
for i=1:length(event)
    val(i) = event(i).value;
    sample(i) = event(i).sample;
end
val1 = unique(val);
Evnt_IDs = mode(val);
%%

if ~isempty(Index)
    cfg                         = [];
    cfg.dataset                 = datafile;
    cfg.trialfun                = 'ft_trialfun_general'; % this is the default
    cfg.trialdef.eventtype      = 'STI101';
    cfg.trialdef.eventvalue     = Evnt_IDs; % the value of the stimulus trigger for fully incongruent (FIC).
    cfg.trialdef.prestim        = 1; % in seconds
    cfg.trialdef.poststim       = 3; % in seconds
    cfg = ft_definetrial(cfg);
    
    cfg.channel = {'MISC001'};
    cfg.hpfreq = 70;
    cfg.demean = 'yes';
    speech_data = ft_preprocessing(cfg);
    
    %% Speech onset detection
    tt=[];
    
    for i=1:length(speech_data.trial)
        tmp = speech_data.trial{1,i} - mean(speech_data.trial{1,i});
        tmp = detrend(tmp);
        tmp1 = zeros(1,length(tmp));
        tmp1(1,200:end-200) = tmp(1,200:end-200);
        
        %     [thres_buf,env, bin] = envelop_hilbert_modified(abs(tmp1),20,1,20,1);
        [thres_buf,env, bin] = envelop_hilbert_modified(abs(tmp1),20,1,20,0);
        
        [a,b] = find(thres_buf > 0.8.*max(thres_buf));
        [d,initCross,finalCross,nextCross,midRef] =  dutycycle(bin);
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
        
        %     figure, plot(speech_data.time{1}, abs(hilbert(speech_data.trial{1,i}))),
        %     hold on,
        %     scaled = (ipoints * (abs(min(speech_data.time{1,i})) + abs(max(speech_data.time{1}))))/length(speech_data.trial{1,i});
        %     scaled = scaled - abs(min(speech_data.time{1,i}));
        %     vline(scaled,'g',['speech onset:', num2str(scaled),'sec']),
        %     box off;
        %     set(gca,'color','none');
        
    end
    
    thre  = speech_data.time{1}(round(0.5.*mean(ipoints_all)));
    if thre < 0.4, thre = 0.4; end
    
    L = length(tt);
    idx_good = find(ipoints_all >= 0.5.*mean(ipoints_all));
    idx_good_1 = find(tt >= 0.4); % responses after 400ms are considered valid.
    idx_good = intersect(idx_good,idx_good_1);
    idx_bad  = find(ipoints_all < 0.5.*mean(ipoints_all));
    
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
    
    
    k=1;
    nfig = ceil(length(speech_data.trial)/12);
    for i=1:nfig
        figure,
        for j=1:12
            if k<=length(speech_data.trial)
                subplot(3,4,j)
                plot(speech_data.time{1}, abs(hilbert(speech_data.trial{1,k}))),
                hold on,
                vline(tt(k),'r',['speech onset:', num2str(tt(k)),'sec']),
                if ~~ipoints_good(k)==1
                    title(num2str(k))
                else
                    title('BAD')
                end
                k=k+1;
            end
        end
    end
    
    stt = std(tt(idx_good));
    mtt = mean(tt(idx_good));
    
    disp([num2str(mtt),'+-', num2str(stt),' Sec'])
    
    
    %%
    savepath = fullfile(outd.sub,'speech');
    if exist(savepath, 'file') == 0, mkdir(savepath), end
    
    print('speech/speech','-dpng');
    badSpeech = table(idx_bad');
    textfile_rej = 'speech/bad_speech';
    badSpeech.Properties.VariableNames{'Var1'} = 'bSpeechs';
    writetable(badSpeech,textfile_rej,'Delimiter',' ');
    
    save(fullfile(savepath,'RT'),'tt','mtt','stt')
    
    %%
    
else
    disp('no speech data available')
end
disp(['Data was saved at,'])
disp(savepath)

%% BS update
%- Input dir
indir = '/MEG_data/epilepsy';
%- Output dir
% outdir = '/MEG_data/Vahab/Processed_data';
bsdatadir = fullfile(indir,subj,'brainstorm_db/data');
cd(bsdatadir)
[subjdir, spath] = uigetfile;

test = load([spath,subjdir]);
test.BadTrials = idx_bad';
save([spath,subjdir],'-struct', 'test'),



%%
% print('ica/ica2','-dpng');


%% BACKUP copies,
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



%%
% clear cln_data
% ic_selection = 1; % 1: manual, 2: automated
% switch task
%     
%     case {1,2}
%         Run_preprocess
%         
%         %- Notch filtering: 30Hz
%         if flag.notch ==1, Run_notch, end %- Notch filtering of 30Hz
%         
%         savetag1 = [tag,'_',subj,'_BS_channels_speech_corrected'];
%         savetag2 = [tag,'_',subj,'_IC_data_speech_corrected'];
%         %- BrainStorm export preprocessed ft-IC
%         
%         %- Speech
%         if flag.speechanalysis ==1
%             Run_speech
%         end
%         
%         %%
%         % export to brainstorm
%         Run_bs
%         
%     case 3
%         Run_preprocess_motor
%         cln_data = clnl_data; % left data condition
%         if flag.notch ==1, Run_notch, end %- Notch filtering of 30Hz
%         savetag1 = [tag,'_',subj,'_left_BS_channels'];
%         savetag2 = [tag,'_',subj,'_left_IC_data'];
%         Run_bs
%         
%         cln_data = clnr_data; % right data condition
%         if flag.notch ==1, Run_notch, end
%         savetag1 = [tag,'_',subj,'_right_BS_channels'];
%         savetag2 = [tag,'_',subj,'_right_IC_data'];
%         Run_bs
% end
% 
% %% Selecting freq of interest
% fmax = 40;
% datain = cln_data;
% if flag.freq == 1
%     stag = 'tsk_baseline';
%     Run_freq
%     disp(['time_of_interest:',num2str(time_of_interest),'sec']);
%     disp(['freq_of_interest:',num2str(freq_of_interest),'Hz']);
%     L = 0.3;
% end
% pause(10)
% 
% %%
% close all
% set(0,'DefaultFigureWindowStyle','normal')
% addpath('/usr/local/MATLAB_Tools/brainstorm3')
% %brainstorm
% 
% %%
% cd (bssavedir)
% clear
% 
% disp('completed, data are ready to import into BS!')




