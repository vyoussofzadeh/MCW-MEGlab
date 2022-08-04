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
%
ft_path = fullfile('/opt/matlab_toolboxes/ft_packages/fieldtrip_20190419');
addpath(ft_path);
ft_defaults

allpath = [];
allpath.ft_path18 = fullfile('/opt/matlab_toolboxes/ft_packages/fieldtrip_041718'); % needed for IC plotting
allpath.ft_path = ft_path;
allpath.script_path = script_path;

% PAC toolbox
path_to_pacMEG      = '/group/prishah/work/tACS/eeg/code/functions';
addpath(genpath(path_to_pacMEG));

%% Load data
taskd = 'tone'; %options: 'pstm' or 'tone'
load(fullfile(ftoutdir, ['all_dataclean_', taskd, '_', subjid, '.mat'])) %loads alldata
ftdatatone = alldata;
clear alldata
taskd = 'pstm'; %options: 'pstm' or 'tone'
load(fullfile(ftoutdir, ['all_dataclean_', taskd, '_', subjid, '.mat'])) %loads alldata
ftdatastm = alldata;

taskd = 'tone'; %options: 'pstm' or 'tone'
load(fullfile(ftoutdir, ['MI', taskd, '_plv_', subjid, '.mat'])) %loads MI
MItone = MI;
clear MI
taskd = 'pstm'; %options: 'pstm' or 'tone'
load(fullfile(ftoutdir, ['MI', taskd, '_plv_', subjid, '.mat'])) %loads MI
MIstm = MI;
clear MI
%% Generate dummy freqanalysis outputs
cfg.output = 'pow';
cfg.method = 'mtmconvol';
cfg.taper = 'hanning';
cfg.foi = amfq;
cfg.t_ftimwin = ones(length(cfg.foi),1).*0.5;
cfg.toi = ftdatastm.time{1}(1:1:7);
cfg.keeptrials = 'yes';
ftfqtone = ft_freqanalysis(cfg, ftdatatone);
ftfqstm = ft_freqanalysis(cfg, ftdatastm);

ftfqtone.powspctrm = atanh(permute(MItone,[2,1,3,4]));
ftfqstm.powspctrm = atanh(permute(MIstm,[2,1,3,4]));

%% Generate neighbor
load('/opt/matlab_toolboxes/ft_packages/fieldtrip_041718/fieldtrip-master/template/neighbours/neuromag306mag_neighb.mat');%neuromag306cmb_neighb.mat')

%% Cluster analysis
cfg = [];
cfg.method = 'montecarlo';
cfg.statistic = 'ft_statfun_indepsamplesT';
cfg.latency = 'all';
cfg.frequency = 'all';
cfg.correctm         = 'cluster';
cfg.clusteralpha     = 0.050;
cfg.clusterstatistic = 'maxsum';
cfg.tail             = 0;
cfg.clustertail      = 0;
cfg.alpha            = 0.025;
cfg.numrandomization = 500;
cfg.neighbours = neighbours; 
cfg.design           = [ones(1,size(ftfqstm.powspctrm,1)), ones(1,size(ftfqtone.powspctrm,1))*2];
cfg.ivar             = 1;

[stat] = ft_freqstatistics(cfg, ftfqstm, ftfqtone);

%% Plot cluster results

close all
%evaluate positive clusters
pos_cluster_pvals = [stat.posclusters(:).prob];
pos_clust = find(pos_cluster_pvals < 0.025);
pos       = ismember(stat.posclusterslabelmat, pos_clust);
%negative clusters
neg_cluster_pvals = [stat.negclusters(:).prob];
neg_clust         = find(neg_cluster_pvals < 0.025);
neg               = ismember(stat.negclusterslabelmat, neg_clust);

allstat = stat.stat;
inclus = pos;
sigclus = inclus.*allstat;
ind = find(sigclus);
[i1,i2,i3]=ind2sub(size(sigclus),ind);
sigch = unique(i1);
%
%alldata.label
%[lid,nn]=ismember('MEG1641',alldata.label)

for ii= 1:length(sigch)
    nn = sigch(ii);
    disp(alldata.label{nn})
    close all
    insig = squeeze(allstat(nn,:,:));
    inmask = squeeze(inclus(nn,:,:));
    figure('color', 'w');
    set(gcf, 'Position', [371 252 874 709])
    pcolor(phfq,amfq,insig); % colormap
    axislim = max(sigclus(:))*0.65;
    caxis([-axislim axislim]) %threshold
    shading interp; colormap(jet);hold on; c =colorbar; %shading, colorbar
    contour(phfq,amfq,inmask,[1],'--','Color','black','LineWidth',3) %stats mask
    set(gca,'FontSize',30);
    ylabel('Amplitude Frequency (Hz)','FontSize',25); xlabel('Phase Frequency (Hz)','FontSize',25) %axis labels
    ylabel(c,'t-value','FontSize',25);
    pause
end

%% Topographic map individual frequencies
stat.raweffect = squeeze(mean(ftfqstm.powspctrm,1)) - squeeze(mean(ftfqtone.powspctrm,1));
[~,p]=ismember(4,phfq);
[~,a]=ismember([24,40],amfq)
fqind1 = p;%phfq(3:6);
fqind2 = a(1):a(2);%[3,8,13];%amfq([3,8,13]);
allstat = stat.stat;
inclus = pos;
sigclus = inclus.*allstat;
ind = find(sigclus);
[i1,i2,i3]=ind2sub(size(sigclus),ind);
sigch = unique(i1);
totch = length(alldata.label);
for k1 = fqind1
    f1 = phfq(k1);
    for k2=fqind2    
        close all
        int = zeros(totch,1);    
        f2 = amfq(k2);
        disp(['plotting for ' num2str(f1) ' and ' num2str(f2)])
        int(1:totch) = inclus(1:totch,k2,k1);       
        data = [];
        data.avg = sigclus(1:totch,k2,k1);
        data.label = stat.label;
        data.time = 1;
        data.dimord = 'chan_time';        
        cfg=[];
        cfg.xlim = 'maxmin';
        cfg.zlim = 'maxmin';        
        cfg.highlight = 'on';
        cfg.highlightchannel = find(int);
        cfg.layout = 'neuromag306mag.lay';        
        cfg.figure= 'gca';        
        ft_topoplotER(cfg,data)
        pause
    end
    
end

%% Topographic map for average within phase and amplitude frequencies
[~,fqind1] = ismember(10,phfq);
phfq(fqind1)
[~,a] = ismember([28,32],amfq);
fqind2 = a(1):a(2);
amfq(fqind2)
allstat = stat.stat;
inclus = pos;
sigclus = inclus.*allstat;
ind = find(sigclus);
[i1,i2,i3]=ind2sub(size(sigclus),ind);
sigch = unique(i1);
totch = length(alldata.label);

close all
int = zeros(totch,1);   
temp = mean(inclus(1:totch,fqind2,fqind1),2);
int(1:totch) = squeeze(mean(temp,3));       
data = [];
temp = mean(sigclus(1:totch,fqind2,fqind1),2);
data.avg = squeeze(mean(temp,3));
data.label = stat.label;
data.time = 1;
data.dimord = 'chan_time';        
cfg=[];
cfg.xlim = 'maxmin';
cfg.zlim = 'maxmin';        
%cfg.highlight = 'on';
%cfg.highlightchannel = find(int);
cfg.layout = 'neuromag306mag.lay';        
cfg.figure= 'gca';        
ft_topoplotER(cfg,data)

%% Analytic analysis

cfg = [];
cfg.method = 'analytic';
cfg.statistic = 'ft_statfun_indepsamplesT';
cfg.latency = 'all';
cfg.frequency = 'all';
cfg.tail             = 0;
cfg.alpha            = 0.001;
cfg.numrandomization = 2000;
cfg.neighbours = neighbours; 
cfg.design           = [ones(1,size(ftfqstm.powspctrm,1)), ones(1,size(ftfqtone.powspctrm,1))*2];
cfg.ivar             = 1;

[stat] = ft_freqstatistics(cfg, ftfqstm, ftfqtone);
%%
pos       = stat.mask;
allstat   = stat.stat;
inclus = pos;
sigclus = inclus.*allstat;
ind = find(sigclus);
[i1,i2,i3]=ind2sub(size(sigclus),ind);
sigch = unique(i1);

for ii= 1:length(sigch)
    nn = sigch(ii);
    disp(alldata.label{nn})
    close all
    insig = squeeze(allstat(nn,:,:));
    inmask = squeeze(inclus(nn,:,:));
    figure('color', 'w');
    set(gcf, 'Position', [371 252 874 709])
    pcolor(phfq,amfq,insig); % colormap
    axislim = max(sigclus(:))*0.65;
    caxis([-axislim axislim]) %threshold
    shading interp; colormap(jet);hold on; c =colorbar; %shading, colorbar
    contour(phfq,amfq,inmask,[1],'--','Color','black','LineWidth',3) %stats mask
    set(gca,'FontSize',30);
    ylabel('Amplitude Frequency (Hz)','FontSize',25); xlabel('Phase Frequency (Hz)','FontSize',25) %axis labels
    ylabel(c,'t-value','FontSize',25);
    pause
end