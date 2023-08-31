%% The Spike Detection MEG pipline

% Spike Detection MEG pipline
% Writtern by MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Lastest update: 08/09/2022

clear; clc, close('all'); warning off

%% FieldTrip toolbox
restoredefaultpath % reset the default path
ft_path ='/opt/matlab_toolboxes/ft_packages/Stable_version/fieldtrip-master';
addpath(ft_path);
ft_defaults

addpath('/data/MEG/Research/awang/Scripts/func')

datadir = '/data/MEG/Research/awang/Epil_annotated_data/annotated_data_anonymized';

% if exist(savedir, 'file') == 0, mkdir(savedir);  end

addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/functions_new')

%%
cd(datadir)
d = rdir([datadir,'/*.mat']);

%%
clear subj run sub_run
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    tkz = tokenize(name,'_');
    subj{i} = [tkz{1}, '_', tkz{2}];
    sub_run{i,:} = [subj{i}, '_', tkz{3}];
end
[sub_run_unq,IA,IC] = unique(sub_run);
% disp(sub_run_unq);

sub_run_unq1 = [];
for i=1:length(sub_run_unq)
    sub_run_unq1{i} = [num2str(i), '_', sub_run_unq{i}];
end
disp(sub_run_unq1')

%%
disp('Enter sub')
subid = input('');

%%

foi = [1,30];

outsum = [];

for i= subid%1:length(d)
    disp([num2str(i),'/',num2str(length(d))])
    [pathstr, name] = fileparts(d(i).name);
    load(d(i).name);
    
    for j=1:length(anot_data_all)
        anot_data = anot_data_all{j};
        
        cfg = [];
        cfg.channel = 'MEG*';
        MEG_data = ft_selectdata(cfg, anot_data);
        cfg.channel = 'EEG*';
        EEG_data = ft_selectdata(cfg, anot_data);
        
        cfg = [];
        cfg.foilim = [1,30];
        cfg.tapsmofrq = 2;
        cfg.taper = 'hanning'; % Tapering window type
        cfg.pad  = 4;
        cfg.plotflag = 1;
        cfg.saveflag = 0;
        [freq, ff, psd,tapsmofrq] = do_fft(cfg, MEG_data); title('meg')
        [freq, ff, psd,tapsmofrq] = do_fft(cfg, EEG_data); title('eeg')
        
        
        %         cfg = [];
        %         cfg.plot = 0;
%         outsum = do_conn(MEG_data.trial{1});
%         figure,imagesc(outsum)
%         figure, plot(mean(outsum))
        
        
        cfg = [];
        cfg.blocksize = anot_data.time{1}(end) - anot_data.time{1}(1);
        cfg.viewmode = 'vertical'; %butterfly';
        cfg.continuous = 'yes';
        cfg.axisfontsize = 7;
        cfg.fontsize = 7;
        cfg.channel = 'EEG*';
        cfg.preproc.demean = 'yes';
        cfg.position = [300   900   500   1500];
        ft_databrowser(cfg, anot_data);
        cfg.channel = 'MEG*';
        cfg.position = [850   900   500   1500];
        ft_databrowser(cfg, anot_data);
        
        pause,
        close all,
    end
end

%%

% cfg = [];
% cfg.toi = [MEG_data.time{1}(end), MEG_data.time{1}(1)];
% cfg.subj = 'test';
% cfg.savefile = [];
% tfr = do_tfr(cfg, MEG_data)
% 
% 
% cfg              = [];
% %         cfg.method       = 'mtmfft';
% cfg.method     = 'mtmconvol';
% cfg.output       = 'pow';
% cfg.foi        = 1:1:10;
% cfg.t_ftimwin  = 3./cfg.foi;
% %         cfg.tapsmofrq  = 0.8 *cfg.foi;
% %         cfg.foilim       = [1,30];
% cfg.toi = 1.9:0.01:2.0;
% %         cfg.tapsmofrq    = 2;
% cfg.taper        = 'hanning';
% cfg.pad          = 4;
% freq2             = ft_freqanalysis(cfg, MEG_data);
% tmp = freq2.powspctrm;
% tmp(isnan (tmp))=0;
% freq2.powspctrm = tmp;
% 
% cfg = []; cfg.savepath = 1; cfg.savefile = [];
% cfg.toi = [MEG_data.time{1}(1), MEG_data.time{1}(end)];
% cfg.bslcorr = 1; cfg.plotflag = 1; cfg.title = 'test';
% cfg.fmax = 30;
% [~,~, tfr_val]    = do_tfr_plot(cfg, freq2);
% 
% cfg = [];
% cfg.foilim = [1,30];
% cfg.tapsmofrq = 2;
% cfg.taper = 'hanning'; % Tapering window type
% cfg.pad  = 4;
% cfg.plotflag = 1;
% cfg.saveflag = 0;
% [freq, ff, psd,tapsmofrq] = do_fft(cfg, MEG_data);
% 
% 
% cfg = [];
% cfg.method = 'mtmconvol'; % Multitaper time-frequency analysis
% cfg.output = 'pow'; % Output power spectra
% cfg.taper = 'hanning'; % Tapering window type
% cfg.foi = 4:20; % Frequency of interest (1 to 30 Hz)
% cfg.toi = 1.9:0.01:2.0; % Time of interest (relative to the segment start)
% cfg.t_ftimwin = 0.2 * ones(size(cfg.foi)); % Length of time window for frequency analysis (in seconds)
% freq = ft_freqanalysis(cfg, MEG_data); % 'data' is your MEG data structure
% 
% cfg = [];
% cfg.colorbar = 'yes'; % Show colorbar
% % cfg.layout = 'your_layout_file.mat'; % Layout file for channel positions
% ft_multiplotTFR(cfg, freq);
% 
% cfg = []; cfg.savepath = 1; cfg.savefile = [];
% cfg.toi = [freq.time(1), freq.time(end)];
% cfg.bslcorr = 1; cfg.plotflag = 1; cfg.title = 'test';
% cfg.fmax = 30;
% [~,~, tfr_val]    = do_tfr_plot(cfg, freq);
% 
% cfg = [];
% cfg.channel = 'MEG*';
% MEG_data = ft_selectdata(cfg, anot_data);
%
%
%         cfg = [];
%         cfg.toi = [MEG_data.time{1}(end), MEG_data.time{1}(1)];
%         cfg.subj = 'test';
%         cfg.savefile = [];
%         tfr = do_tfr(cfg, anot_data)
%
%
%
%         cfg = [];
%         cfg.output     = 'pow';
%         cfg.channel    = 'all';
%         cfg.method     = 'mtmconvol';
%         cfg.method     = 'wavelet';
%         if foi(2) - foi(1) < 10
%             cfg.foi        = foi(1):1:foi(2);
%         else
%             cfg.foi        = foi(1):2:foi(2);
%         end
%         cfg.keeptrials = 'yes';
%         cfg.t_ftimwin  = 3 ./ cfg.foi;
%         cfg.tapsmofrq  = 0.8 * cfg.foi;
%         cfg.toi        = anot_data.time{1}(1):0.05:anot_data.time{1}(end);
%         tfr_data        = ft_freqanalysis(cfg, anot_data);
% %         tfr_data.powspctrm = squeeze(tfr_data.powspctrm);
%
%         cfg = []; cfg.savepath = 1; cfg.savefile = [];
%         cfg.toi = [tfr_data.time(1), tfr_data.time(end)];
%         cfg.bslcorr = 1; cfg.plotflag = 1; cfg.title = 'test';
%         cfg.fmax = 40;
%         [~,~, tfr_val]    = do_tfr_plot(cfg, tfr);