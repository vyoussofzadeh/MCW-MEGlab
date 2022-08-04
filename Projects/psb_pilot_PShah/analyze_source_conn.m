%% The psb_pilot

% MEG phase amplitude coupling
% Writtern by MCW group, Shah-Basak, Priyanka <prishah@mcw.edu>
% input data from MEG (pre-) processing pipeline by Youssofzadeh, Vahab
% Lastest update: 05/25/2022

clear; clc, close('all'); warning off

%% Paths
restoredefaultpath
maindir = '/group/prishah/LanguageMEG/psb_pilot';
ssid = 'ss_pilot_2';
subjid = 'pilot2';
datcol = '220512';

script_path = maindir;%'/data/MEG/Research/psb_pilot';
addpath(genpath(script_path));
%- Input dir
indir = fullfile(maindir, ssid, datcol, 'tsss');%'/data/MEG/Research/psb_pilot/ss_pilot_2/220512/tsss';
inpath = fullfile(indir, 'STM_Block1_raw.fif');

%- Output dir
%outdir = fullfile(maindir, 'FT');
ftoutdir = fullfile(maindir, ssid, datcol, 'meg');
%- MRI dir
mridir = fullfile(maindir, ssid, datcol, 'mri');
ftanatdir = fullfile(maindir, ssid, datcol, 'meg' , 'anat');
%
ft_path = fullfile('/opt/matlab_toolboxes/ft_packages/fieldtrip-20210517'); %fieldtrip_20190419
addpath(ft_path);
ft_defaults
%addpath('/opt/matlab_toolboxes/ft_packages/fieldtrip-20210517/')
%ft_defaults
allpath = [];
allpath.ft_path18 = fullfile('/opt/matlab_toolboxes/ft_packages/fieldtrip_041718'); % needed for IC plotting
allpath.ft_path = ft_path;
allpath.script_path = script_path;

%- add BrainNet Viewer and BCT path
addpath(genpath('/group/prishah/work/tACS/toolboxes/BrainNetViewer'))
bnetpath = '/group/prishah/LanguageMEG/bnet';

addpath('/group/prishah/work/tACS/toolboxes/2019_03_03_BCT')
%% Acquire aal atlas information
disp('acquire aal atlas info')
aal90 = readtable('/group/prishah/LanguageMEG/bnet/aal90.txt', 'FileType', 'text', 'ReadVariableNames', false,'HeaderLines', 0);
aal90.Properties.VariableNames = { 'x', 'y', 'z','lobe','values', 'labels'};  
nn = size(aal90,1);
% aal90(ia{ss}(ia{ss}<=90),:) = [];
l2l = 1:2:nn;
r2r = 2:2:nn;
aal90 = [aal90(l2l,:); aal90(r2r,:)];
aal90_coords = aal90(:,1:4);
aal90_labels = aal90(:,6);

%% Load connectivity matrices 
close all;
fq  ='8'; %frequency of interest options
taskd = 'tone'; %options: 'pstm' or 'tone' 
load(fullfile(ftoutdir, ['imagcoh_', taskd, '_', fq ,'Hz_aal.mat'])) %loads wplitone
conntone=[]; conntone = abs(icohtone.cohspctrm);

taskd = 'pstm'; %options: 'pstm' or 'tone'
load(fullfile(ftoutdir, ['imagcoh_', taskd, '_', fq ,'Hz_aal.mat'])) %loads wplistm
connstm=[]; connstm = abs(icohstm.cohspctrm);

%% Rearrange connectivity matrices with left nodes first and then right nodes
nlabel =wplitone.label;
l2l = find(contains(nlabel,'_L'));
r2r = find(contains(nlabel,'_R'));
nlabel = [nlabel(l2l); nlabel(r2r)];    
pr1 = find(contains(nlabel, 'Paracentral_Lobule_R'));
nlabel(pr1(1))=[];
l2l(pr1(1)) = [];
tconnstm = [connstm(l2l, l2l), connstm(l2l, r2r); connstm(r2r, l2l), connstm(r2r, r2r)];
tconntone = [conntone(l2l, l2l), conntone(l2l, r2r); conntone(r2r, l2l), conntone(r2r, r2r)];

avgDiff = tconnstm-tconntone;
avgDiff(isnan(avgDiff))=888;
avgDiff=avgDiff.*~(eye(size(avgDiff)));
%% Visualize using brainnet viewer
top = 0.70;%
bavgDiff = avgDiff;
bavgDiff(bavgDiff<0)=0;
bavgDiff = abs(bavgDiff);
bavgDiff(bavgDiff<=max(max(bavgDiff))*top) = 0;
figure('color','white')
subplot(121); imagesc(avgDiff); colorbar; yline(32); xline(32); daspect([1 1 1]); colormap('jet'); 
subplot(122); imagesc(bavgDiff); colormap('jet'); colorbar; daspect([1 1 1]); title(['top ' num2str((1-top)*100) '%'])
pause;
path_base = [bnetpath '/avg_conn']; 
dlmwrite([path_base '-up.edge'], bavgDiff, 'delimiter', ' ');

disp('acquire aal atlas info')
aal90 = readtable([bnetpath '/aal90.txt'], 'FileType', 'text', 'ReadVariableNames', false,'HeaderLines', 0);
aal90.Properties.VariableNames = { 'x', 'y', 'z','lobe','values', 'labels'};  
% aal90(ia{ss}(ia{ss}<=90),:) = [];
l2l = 1:2:nn;
r2r = 2:2:nn;
aal90 = [aal90(l2l,:); aal90(r2r,:)];
aal90_coords = aal90(:,1:4);
aal90_labels = aal90(:,6);

ns = strengths_und(bavgDiff);
tab_node_up = [aal90_coords, array2table(ns'), aal90_labels];
writetable(tab_node_up, [path_base '-up.node'],...
    'WriteVariableNames', false, 'FileType', 'text', 'Delimiter', '\t');
BrainNet_MapCfg(['/group/prishah/work/tACS/toolboxes/BrainNetViewer/Data/SurfTemplate/BrainMesh_ICBM152.nv'], ...
                [path_base '-up.edge'], ...
                [path_base '-up.node'], ... %'-eigcent-up.node'                        
                [bnetpath '/config.mat']);
            

