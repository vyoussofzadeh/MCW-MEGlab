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

datadir = '/data/MEG/Research/SpikeDectection/Epil_annotated_data/annotated_data_anonymized';

% if exist(savedir, 'file') == 0, mkdir(savedir);  end

addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/helper')

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
for i= subid%1:length(d)
    disp([num2str(i),'/',num2str(length(d))])
    [pathstr, name] = fileparts(d(i).name);
    load(d(i).name);
    
    for j=1:length(anot_data_all)
        anot_data = anot_data_all{j};
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