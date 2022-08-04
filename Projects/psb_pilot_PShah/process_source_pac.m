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
indir = fullfile(maindir, ssid, datcol, 'tsss');%for sensor-level data

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

% PAC toolbox
path_to_pacMEG      = '/group/prishah/work/tACS/eeg/code/functions';
addpath(genpath(path_to_pacMEG));

%% Load data
taskd = 'tone'; %options: 'pstm' or 'tone'
load(fullfile(ftoutdir, ['source_alltrials_' taskd '.mat'])) %loads output
alldata = virttone; %csourcestm;
clear virttone %csourcestm
fs = 2000;
nchannel = length(alldata);
insrc=find(~cellfun(@isempty,alldata));
ntrial = size(alldata{insrc(1)},1);

%% Run PAC analysis
MI = zeros(ntrial, nchannel, length(amfq), length(phfq));

for chn = 1:nchannel
    
    if ismember(chn,insrc)
        disp(['Computing channel: ' num2str(chn)])
        %PAC specs
        cfg                     = [];
        cfg.Fs                  = fs;
        cfg.phase_freqs         = phfq;
        cfg.amp_freqs           = amfq;
        cfg.method              = 'canolty';%list_of_methods{method};
        cfg.filt_order          = 3;
        %cfg.mask               = [691 1051];
        %cfg.surr_method         = 'swap_blocks';
        cfg.surr_N              = 200;
        cfg.amp_bandw_method    = 'number';
        cfg.amp_bandw           = 15;
        %real data

        for tt = 1:ntrial
            PAC_signal=[]; PAC_signal = alldata{chn}(tt,:);
            %PAC_signal=[]; PAC_signal = squeeze(alldata.trial(tt,chn,:))';                
            MI(tt,chn,:,:)          = PACmeg(cfg,PAC_signal);
        end%tt trial
    end
    %mMIreal{ss}(chn,:,:) = squeeze(mean(MIreal,1));
   
end %chn channels
save(fullfile(ftoutdir,['MI' taskd '_source_' subjid '.mat']),'MI')

%% Plotting PAC
[lid,nn]=ismember('SupraMarginal_L',alldata.label)
if lid    
    indata = squeeze(mean(MI(:,nn,:,:),1));
    plot_comod(phfq,amfq,indata); colorbar;
end


figure('color', 'w');
pcolor(phfq,amfq,indata); % colormap
axislim = max(indata(:))*0.65;
caxis([-axislim axislim]) %threshold
shading interp; colormap(jet);hold on; c =colorbar; %shading, colorbar
%contour(phfq,amfq,posclus,[1],'--','Color','black','LineWidth',3) %stats mask
set(gca,'FontSize',30);
ylabel('Amplitude Frequency (Hz)','FontSize',25); xlabel('Phase Frequency (Hz)','FontSize',25) %axis labels
ylabel(c,'t-value','FontSize',25);