%% The Spike Detection MEG pipline

% Spike Detection MEG pipline
% Writtern by MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Lastest update: 08/25/2022

clear; clc, close('all'); warning off

%% FieldTrip toolbox
restoredefaultpath % reset the default path
% ft_path ='/opt/matlab_toolboxes/ft_packages/Stable_version/fieldtrip-master';
ft_path = '/opt/matlab_toolboxes/ft_packages/fieldtrip_latest';
addpath(ft_path);
ft_defaults

addpath('/data/MEG/Research/awang/Scripts/func')

datadir = '/data/MEG/Research/SpikeDectection/Epil_annotated_data/raw_data';

addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/helper')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/Deeplearning_spike/Squiggles/func')

addpath(genpath('/opt/matlab_toolboxes/mne_matlab/matlab'))

% if exist(savedir, 'file') == 0, mkdir(savedir);  end

%%
cd(datadir)
d = rdir([datadir,'/**/*.fif']);

%%
clear subj run sub_run datafile
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    tkz = tokenize(pathstr,'_'); 
    tkz1 = tokenize(tkz{end-2},'/');
    sub_run{i} = [num2str(i), ':', tkz1{2}, '_', tkz{end-1},'_', tkz{end}];
    datafile{i} = d(i).name;
end
% [sub_run_unq,IA,IC] = unique(sub_run);
% disp(sub_run_unq')
disp(sub_run')

%%
disp('Enter data sub')
subid = input('');

%%
disp('Select frequncy range (Hz)');
disp('1: Wideband 4-40');
disp('2: Slow-rate 1-5')
disp('3: High-rate 20-40');
disp('4: Other frequncy');
ask.freq_occur_sel = input(':');

switch ask.freq_occur_sel
    case 1
        foi = [4,40];
    case 2
        foi = [1,5];
    case 3
        foi = [25,40];
    case 4
        disp('enter range [f1,f2]Hz:'); freq_range_sel = input(':');
        foi = [freq_range_sel(1), freq_range_sel(2)];
end

%%
cfg = []; cfg.layout = 'neuromag306mag.lay'; lay = ft_prepare_layout(cfg);

for i= subid %1:length(d)
    
    disp([num2str(i),'/',num2str(length(d))])
    datafile = d(i).name;
    
    %- Preprocess data
    cfg = []; cfg.datafile = datafile;
    cfg.modal = 'meg';
    [ds_meg, aft_meg] = do_preprocess_spike_data (cfg, datafile);
    cfg.modal = 'eeg';
    [ds_eeg, aft_eeg] = do_preprocess_spike_data (cfg, datafile);
    
    %%
    cfg = [];
    cfg.channel = 'MEG*';
    sel_meg = ft_selectdata(cfg,ds_meg);
    cfg.channel = 'EEG*';
    sel_eeg = ft_selectdata(cfg,ds_eeg);
    
    cfg          = [];
    cfg.hpfilter = 'yes';
    cfg.lpfilter = 'yes';
    cfg.hpfiltord = 3;
    cfg.hpfreq = foi(1);
    cfg.lpfreq = foi(2);
    f_megdata = ft_preprocessing(cfg, sel_meg);
    f_eegdata = ft_preprocessing(cfg, sel_eeg);
    
    cfg = [];
    cfg.aftval = 0;
    cfg.aft = aft_eeg;
    af_eegdata = do_appyaft(cfg,f_eegdata);
    cfg.aft = aft_meg;
    af_megdata = do_appyaft(cfg,f_megdata);
        
    cfg = [];
    cfg.blocksize = 10;
    cfg.viewmode = 'vertical'; %butterfly';
    cfg.continuous = 'yes';
    cfg.axisfontsize = 7;
    cfg.fontsize = 7;
    cfg.channel = 'EEG*';
    %     cfg.preproc.demean = 'yes';
    cfg.position = [300   900   500   1500];
    ft_databrowser(cfg, af_eegdata);
    cfg.channel = 'MEG*';
    cfg.position = [850   900   500   1500];
    ft_databrowser(cfg, af_megdata);
    
end

%%
pause,
sdir = '/data/MEG/Research/awang/Epil_annotated_data/sample_testdata';
close all
disp('save data (y: yes)?')
asksavedata = input('','s');
switch asksavedata
    case 'y'
        save(fullfile(sdir,[sub_run_unq{subid},'.mat']),'af_megdata','af_eegdata')
end
cd(sdir)

%%
clear

