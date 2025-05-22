%% The Spike Detection MEG pipline

% Spike detection MEG pipline
% Writtern by MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Lastest update: 08/09/2022

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
flag.analysis = 1;

%% Datalog (subject details)
Datalog = [];

%% Initial settings
cd '/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/FT';
restoredefaultpath
cd_org = cd;
addpath(genpath(cd_org));

%- Input dir
indir = '/MEG_data/Research_studies/Epil_annotated_data/annotated_info';
%- Output dir
outdir = '/MEG_data/Research_studies/Epil_annotated_data/annotated_data';

%- Adding path
cfg_init = [];
cfg_init.path_tools = '/usr/local/MATLAB_Tools';
[allpath, atlas] = vy_init(cfg_init);

%%
cd(indir)
d = rdir([indir,'/*.mat']);

%%
clear subj run sub_run
for i=1:length(d)
    [pathstr, name] = fileparts(d(i).name);
    tkz = tokenize(name,'_');
    subj{i} = [tkz{1}, '_', tkz{2}];
    sub_run{i,:} = [subj{i}, '_', tkz{3}];
end
[sub_run_unq,IA,IC] = unique(sub_run);
disp(sub_run_unq);

[sub_unq,IA_s,IC_s] = unique(subj);
disp(sub_unq);

%% Read, preprocess, and save data
if exist(outdir, 'file') == 0, mkdir(outdir), end

for i=1:length(d)
    disp([num2str(i),'/',num2str(length(d))])
    [pathstr, name] = fileparts(d(i).name);
    load(d(i).name);
    if ~exist([fullfile(outdir,name), '.mat'], 'file') && exist(Anot.filename, 'file')
        %     datafile = Anot.filename;
        disp('preprocessing ...')
        cfg = [];
        cfg.dataset = Anot.filename;
        %     cfg.channel = {'megmag', 'meggrad', 'eeg','eog','ecg'};
        cfg.channel = {'megmag', 'meggrad', 'eeg'};
        raw_data = ft_preprocessing(cfg);
        
        cfg = [];
        cfg.resamplefs = 500;
        dsample_data = ft_resampledata(cfg, raw_data);
        
        anot_data_all = []; k=1;
        for j=1:size(Anot.T_tint,1)
            cfg = [];
            if  Anot.T_tint(j,2) -  Anot.T_tint(j,1) > 0.2
                cfg.toilim = Anot.T_tint(j,:)./(raw_data.fsample/dsample_data.fsample);
                anot_data = ft_selectdata(cfg,dsample_data);
                anot_data_all(k,:,:) = anot_data.trial{1}(:,1:67);
                k=k+1;
            end
        end
        
        D = [];
        D.hdr = anot_data.hdr;
        D.label = anot_data.label;
        D.trial = anot_data_all;
        D.Anot = Anot;
        D.time = anot_data.time{1}(1:67);
               
        tkz = tokenize(name,'_');
        save(fullfile(outdir,name), 'D')
%         save(fullfile(outdir,name), 'D', '-v7.3')
    end
end

%% Save orignal MEG data
outdir_rawMEG = '/MEG_data/Research_studies/Epil_annotated_data/raw_data';

if exist(outdir_rawMEG, 'file') == 0, mkdir(outdir_rawMEG), end

flag.plot = 0;

for i=1:length(d)
    disp([num2str(i),'/',num2str(length(d))])
    [pathstr, name] = fileparts(d(i).name);
    load(d(i).name);
    [pathstr2, name2] = fileparts(Anot.filename);
    if ~exist(fullfile(outdir_rawMEG,name,[name2,'.fif']), 'file')
        if exist(fullfile(outdir_rawMEG,name), 'file') == 0, mkdir(fullfile(outdir_rawMEG,name)), end
        copyfile(Anot.filename, fullfile(outdir_rawMEG,name,[name2,'.fif']));
        hdr = ft_read_header(Anot.filename);
        fs = hdr.orig.sfreq;
        first_samp = double(hdr.orig.raw.first_samp);
        
        savefile  = 'evt.txt';
        textfile = fullfile(outdir_rawMEG,name,savefile);
        fid=fopen(textfile,'w');
        event = [first_samp, first_samp/fs, 0, 0];
        fprintf(fid, '%d\t', event); % marker value
        fprintf(fid,'%6s %12s\n','test');
        for j=1:length(Anot.T_tint)
            event(j+1,1) = round(Anot.T_tint(j)*fs) + event(1,1);
            event(j+1,2) = Anot.T_tint(j) + event(1,2);
            event(j+1,3) = 0;
            event(j+1,4) = 5555;
            fprintf(fid,'\n');
            fprintf(fid, '%d\t', event(j+1,:)); % marker value
            fprintf(fid,'%6s %12s\n', 'spk');
        end
        fclose(fid);true
    end
    
    if flag.plot == 1
        %- mne-browse
        cd(fullfile(outdir_rawMEG,name))
        command = ['mbrowse ',[name2,'.fif']];
        system(command);
    end
end

%% Copy to Squiggles
% cd ('/MEG_data/Research_studies/Epil_annotated_data')
% command = 'scp -r /MEG_data/Research_studies/Epil_annotated_data vyoussofzadeh@squiggles.rcc.mcw.edu:/data/MEG/Research/awang';
% system(command)