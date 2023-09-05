%% The Spike Detection MEG pipline

% Spike Detection MEG pipline
% Writtern by MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Lastest update: 08/31/2023


clear; clc, close('all'); warning off

%% FieldTrip toolbox
ft_path ='/opt/matlab_toolboxes/ft_packages/latest/fieldtrip-master';
addpath(ft_path);
ft_defaults

addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/helper')

ddir = '/data/MEG/Research/SpikeDectection/Epil_annotated_data';

%%
disp('1) Spike data')
disp('2) no-Spike data')
datasel = input('');

switch datasel
    case 1
        datadir = fullfile(ddir,'annotated_data');
        savedir = fullfile(ddir,'annotated_data_anonymized');
    case 2
        datadir = fullfile(ddir,'annotated_data_nospike');
        savedir = fullfile(ddir,'annotated_data_nospike_anonymized');
end
if exist(savedir, 'file') == 0, mkdir(savedir);  end

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
disp(sub_run_unq);

%%
cd(savedir)
for i=1:length(d)
    disp([num2str(i),'/',num2str(length(d))])
    [pathstr, name] = fileparts(d(i).name);
    tkz = tokenize(name,'_');
    newName = [tkz{1}(1),'_',tkz{2}(1:2),'_',tkz{3}];
    copyfile(fullfile(datadir,[name,'.mat']), fullfile(savedir,[newName, '.mat']));
end
