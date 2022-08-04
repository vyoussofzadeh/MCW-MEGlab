%% The psb_pilot

% MEG phase amplitude coupling
% Writtern by MCW group, Shah-Basak, Priyanka <prishah@mcw.edu>
% input data from MEG (pre-) processing pipeline by Youssofzadeh, Vahab
% Lastest update: 05/25/2022

clear; clc, close('all'); warning off
%% PAC variables
% freq and time info
phfq = 2:2:14;%
amfq = 16:2:48;

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
%- Output dir
outdir = fullfile(maindir, 'FT');
ftoutdir = fullfile(maindir, ssid, datcol, 'meg');
%- MRI dir
mridir = fullfile(maindir, ssid, datcol, 'mri');
ftanatdir = fullfile(maindir, ssid, datcol, 'meg' , 'anat');
%
ft_path = fullfile('/opt/matlab_toolboxes/ft_packages/fieldtrip-20210517'); %fieldtrip_20190419
addpath(ft_path);
ft_defaults

allpath = [];
allpath.ft_path18 = fullfile('/opt/matlab_toolboxes/ft_packages/fieldtrip_041718'); % needed for IC plotting
allpath.ft_path = ft_path;
allpath.script_path = script_path;

%- PAC toolbox
path_to_pacMEG      = '/group/prishah/work/tACS/eeg/code/functions';
addpath(genpath(path_to_pacMEG));

%- add BrainNet Viewer and BCT path
addpath(genpath('/group/prishah/work/tACS/toolboxes/BrainNetViewer'))
bnetpath = '/group/prishah/LanguageMEG/bnet';
%% Load data

load(fullfile(ftoutdir, 'sample_source_singletrial.mat')) %sourcestm

taskd = 'tone'; %options: 'pstm' or 'tone'
load(fullfile(ftoutdir, ['MI', taskd, '_source_', subjid, '.mat'])) %loads MI
MItone = MI;
clear MI

taskd = 'pstm'; %options: 'pstm' or 'tone'
load(fullfile(ftoutdir, ['MI', taskd, '_source_', subjid, '.mat'])) %loads MI
MIstm = MI;
clear MI

load(fullfile(ftanatdir, 'mri_resliced.mat'))
%% dummy source variables

pval = 12;
aval = [18,28];
[~,p]=ismember(pval,phfq);
[~,a]=ismember(aval,amfq);
disp(['extracting values for phase = ' num2str(pval) ' and amplitude = ' num2str(aval(1)) '-' num2str(aval(2)) ' Hz'])
srcstm = sourcestm;
if size(MIstm,1) == length(sourcestm.trial)    
    for tt = 1:size(MIstm,1)        
        indata=[]; indata = (squeeze(mean(MIstm(tt,:,a(1):a(2),p),3)))';
        srcstm.pow(tt,:) = indata;
    end
end
srcstm=rmfield(srcstm,'trial');

srctone = sourcestm;
if size(MItone,1)<length(sourcestm.trial) 
    ntt = size(MItone,1);
    srctone.cumtapcnt = srctone.cumtapcnt(1:ntt);        
    srctone.df = ntt;
    srctone.trialinfo = srctone.trialinfo(1:ntt);
end
for tt = 1:size(MItone,1)
        indata=[]; indata = squeeze(mean(MItone(tt,:,a(1):a(2),p),3))';
        srctone.pow(tt,:) = indata;
end
srctone = rmfield(srctone,'trial');

srcstm.pow = atanh(srcstm.pow);
srctone.pow = atanh(srctone.pow);
%% run analysis
cfg = [];
cfg.dim = srcstm.dim;
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_indepsamplesT';
cfg.parameter = 'pow';
cfg.correctm         = 'cluster';
cfg.tail             = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
cfg.design(1,:)           = [1:length(srcstm.trialinfo) 1:length(srctone.trialinfo)];
cfg.design(2,:)           = [ones(1,length(srcstm.trialinfo)), ones(1,length(srctone.trialinfo))*2];
%cfg.uvar             = 1;
cfg.ivar             = 2;
[stat] = ft_sourcestatistics(cfg, srcstm,srctone);
%%
%evaluate positive clusters
pos_cluster_pvals = [stat.posclusters(:).prob];
pos_clust = find(pos_cluster_pvals < 0.025);
pos       = ismember(stat.posclusterslabelmat, pos_clust);
stat.pos_clust = pos.*stat.stat;
%negative clusters
neg_cluster_pvals = [stat.negclusters(:).prob];
neg_clust         = find(neg_cluster_pvals < 0.025);
neg               = ismember(stat.negclusterslabelmat, neg_clust);

avgtone = squeeze(mean(MItone,1));
avgstm = squeeze(mean(MIstm,1));
avgdiff = avgstm-avgtone;
stat.avgdiff = squeeze(mean(avgdiff(:,a(1):a(2),p),2));
stat
%% interpolate
cfg = [];
%cfg.downsample = 2;
cfg.parameter= 'pos_clust';
statplot = ft_sourceinterpolate(cfg,stat,mri_resliced);
%%
outfname = fullfile(ftanatdir, ['stat_pac_', num2str(pval), '_' num2str(aval(1)) '-' num2str(aval(2)) 'Hz_STMvsTone']);
[fdir, fname] = fileparts(outfname)
cfg =[];
cfg.filename = outfname;
cfg.filetype = 'nifti';
cfg.parameter = 'pos_clust';
ft_sourcewrite(cfg,ft_convert_units(statplot,'mm'))

system(['3dcalc -prefix ' fdir '/mask_' fname '.nii.gz -a ' fdir '/brain_mask.nii -b ' outfname '.nii -expr ' char(39) '(a*b)' char(39)])
system(['rm ' outfname '.nii'])
%% Average across trials
aalvd = false;
if aalvd
    [lid,nn]=ismember('Angular_L',ftdatastm.label)    
    if lid    
        indata = squeeze(avgdiff(nn,:,:));
        %plot_comod(phfq,amfq,indata); colorbar;
    end
end

indata = squeeze(avgdiff(200,:,:));
figure('color', 'w');
pcolor(phfq,amfq,indata); % colormap
axislim = max(indata(:))*0.65;
caxis([-axislim axislim]) %threshold
shading interp; colormap(jet);hold on; c =colorbar; %shading, colorbar
%contour(phfq,amfq,posclus,[1],'--','Color','black','LineWidth',3) %stats mask
set(gca,'FontSize',30);
ylabel('Amplitude Frequency (Hz)','FontSize',25); xlabel('Phase Frequency (Hz)','FontSize',25) %axis labels
ylabel(c,'MI','FontSize',25);
%%


%% extract MI and  visualize with brainnet
af = find(ismember(amfq,24))
pf = find(ismember(phfq,6))
tone = squeeze(avgtone(:,af,pf));
pstm = squeeze(avgstm(:,af,pf));
diffmi = pstm-tone;

disp('acquire aal atlas info')
aal90 = readtable([bnetpath '/aal90.txt'], 'FileType', 'text', 'ReadVariableNames', false,'HeaderLines', 0);
aal90.Properties.VariableNames = { 'x', 'y', 'z','lobe','values', 'labels'};  
% aal90(ia{ss}(ia{ss}<=90),:) = [];
nn = size(aal90,1);
l2l = 1:2:nn;
r2r = 2:2:nn;
aal90 = [aal90(l2l,:); aal90(r2r,:)];
aal90_coords = aal90(:,1:4);
aal90_labels = aal90(:,6);
diffmi = [diffmi(l2l,:); diffmi(r2r,:)];
tab_node_up = [aal90_coords, array2table(diffmi), aal90_labels];
path_base = [bnetpath '/mi']; 
writetable(tab_node_up, [path_base '-up.node'],...
    'WriteVariableNames', false, 'FileType', 'text', 'Delimiter', '\t');

BrainNet_MapCfg(['/group/prishah/work/tACS/toolboxes/BrainNetViewer/Data/SurfTemplate/BrainMesh_ICBM152.nv'], ...                
                [path_base '-up.node'], ... %'-eigcent-up.node'                        
                [bnetpath '/config.mat']);
