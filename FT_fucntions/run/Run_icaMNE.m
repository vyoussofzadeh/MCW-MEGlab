%% ICA analysis
% pause, close all,
foi = input('enter foi hz, e.g., [18,23]:');

% remove the DC-component
cfg        = [];
cfg.demean = 'yes';
cfg.dftfilter = 'yes';
cfg.hpfilter = 'yes';
cfg.lpfilter = 'yes';
cfg.hpfiltord = 3;
cfg.hpfreq = foi(1);
cfg.lpfreq = foi(2);
%         cfg.channel = {'megmag', 'meggrad'};
fcln_data        = ft_preprocessing(cfg, datain);

%%
cfg            = [];
cfg.method     = 'runica';
cfg.numcomponent = 20;       % specify the component(s) that should be plotted
% cfg.numcomponent = 1;       % specify the component(s) that should be plotted
comp           = ft_componentanalysis(cfg, fcln_data);

%         cfg            = [];
%         cfg.method     = 'pca';
%         cfg.numcomponent = 10;       % specify the component(s) that should be plotted
%         % cfg.numcomponent = 1;       % specify the component(s) that should be plotted
%         comp           = ft_componentanalysis(cfg, cln_data);

%%
%         restoredefaultpath
%         addpath(genpath(allpath.ft18));
%         addpath(genpath(allpath.hcp_path));
%         addpath(genpath(allpath.cd_org));

% cfg = [];
% cfg.viewmode = 'component';
% cfg.layout = lay;
% ft_databrowser(cfg, comp);
% set(gcf, 'Position', [600   400   1200   800]);
% colormap(brewermap(256, '*RdYlBu'));
% title(subj)
% %         set(gcf,'name',subj,'numbertitle','off')

%%
outputsourcedir = fullfile(outd.sub,'mne_ica'); % output dir
if exist(outputsourcedir, 'file') == 0, mkdir(outputsourcedir); end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare the component data in order for ft_sourceanalysis to be able to
% swallow it
mixing   = comp.topo;
channels = comp.topolabel;
% normalisation of the topographies
for i = 1:size(mixing, 2)
    val(i) = 0.01*max(abs(mixing(:, i)));
    mixing(:, i) = mixing(:, i)/val(i);
end

% create a 'timelock' structure
tlck = [];
tlck.label = channels;
tlck.cov = eye(numel(tlck.label));
tlck.time=1;
tlck.grad = datain.grad;
tlck.dimord = 'chan_time';

%%
save(fullfile(outputsourcedir, ['ICA_mixing_', num2str(foi(1)), '_', num2str(foi(2)), 'Hz_', subj, '.mat']), 'tlck','mixing');

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the forward solution

cfg = [];
cfg.vol = individual_headmodel;
cfg.grid = individual_grid;
cfg.grad = cln_data.grad;
cfg.channel = channels;
cfg.normalize = 'yes';
cfg.reducerank = 2;
gridLF = ft_prepare_leadfield(cfg);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do an MNE with different regularisation for each component

% this parameter is hard-coded here
noise_level = 8;

% specify the static part of the cfg for the source reconstruction
cfg               = [];
cfg.method        = 'mne';
cfg.grid          = gridLF;
cfg.vol           = individual_headmodel;

%         cfg.grid = individual_grid;
%         cfg.headmodel = individual_headmodel;
cfg.channel       = channels;
cfg.mne.prewhiten = 'yes';
cfg.mne.noisecov  = eye(numel(channels))*noise_level;

% loop over components, due to component-specific regularisation
for i=1:size(mixing,2)
    
    % use the channel-level topography of the current component
    tlck.avg = mixing(:,i);
    
    % estimate the snr of the current component
    cfg.mne.snr = sqrt(mean((mixing(:,i)-mean(mixing(:,i))).^2))/noise_level;
    noisevec(i) = cfg.mne.snr;
    
    tmp = ft_sourceanalysis(cfg, tlck);
    inside_indices = find(tmp.inside(:))';
    if i==1
        % create the output source structure in the first iteration
        source=tmp;
    else
        % concatenate the reconstructed source level topography to the previously computed ones
        for k = 1:numel(inside_indices)
            source.avg.mom{inside_indices(k)} = cat(2,source.avg.mom{inside_indices(k)}, tmp.avg.mom{inside_indices(k)});
        end
        source.avg.pow = horzcat(source.avg.pow,tmp.avg.pow);
    end
end
%%
clc
%         figure,bar(noisevec), title('SNR')
%         disp(find(noisevec == max(noisevec)))

%%
% pause,
close all,
cfg = [];
cfg.viewmode = 'component';
cfg.layout = lay;
ft_databrowser(cfg, comp);
set(gcf, 'Position', [200   400   1200   800]);
colormap(brewermap(256, '*RdYlBu'));
title(subj)

%%
nIC = input('enter ICs:');
% nIC = 1:5;

for i=1:length(nIC)
    
    cd(outputsourcedir)
    saveID = [subj, '_IC:', num2str(nIC(i))];
    
    source1 = [];
    source1.pow = source.avg.pow(:,nIC(i));
    
    %% subject mapping
    %     source1.dim    = gridLF.dim;
    %     source1.pos = gridLF.pos;
    %     source1.inside  = gridLF.inside;
    %
    %     cfg = [];
    %     cfg.mask = 'pow';
    %     cfg.loc = 'max';
    %     cfg.template = mri_realigned;
    %     cfg.savefile = [];
    %     cfg.volnorm     = 2; % yes: 1
    %     cfg.nslices = 16;
    %     cfg.method  = 'ortho';
    % %     cfg.method  = 'slice';
    %     vy_source_plot(cfg, source1);
    %     set(gcf,'name',['IC:', num2str(nIC(i))],'numbertitle','off')
    %     print([saveID,'_slc'],'-dpng');
    
    %% template/surface mapping
    s_in = source1; mask = 'pow'; thre = 0.3;
    %     Run_surfacemap_template
    
    s_in.pos     = template_grid.pos;
    s_in.dim     = template_grid.dim;
    s_in.inside  = template_grid.inside;
    
    cfg = [];
    cfg.mask = mask;
    cfg.loc = 'max';
    cfg.template = template_mri;
    %     cfg.savefile = [saveID,'_tmp'];
    cfg.savefile = [];
    cfg.volnorm  = 2; % yes: 1
    cfg.nslices = 16;
    cfg.method        = 'ortho';
    source_int = vy_source_plot(cfg, s_in);
    set(gcf,'name',['IC:', num2str(nIC(i))],'numbertitle','off')
    
    %     sIC = input('save surfaces?');
    %     cfg = [];
    %     %     cfg.subj = [];
    %     cfg.subj = saveID;
    %     cfg.mask = mask;
    %     cfg.thre = thre;
    %     if sIC==1
    %         cfg.saveflag = 1;
    %     else
    %         cfg.saveflag = 2;
    %     end
    %     cfg.colorbar = 2;
    %     vy_mapvisualisation(cfg, source_int);
    
end

%% surface mapping
% nIC = input('enter ICs for surface mapping:');
nIC = 1:5;
% sIC = input('save surfaces?');
sIC = 2;
close all,
for i=1:length(nIC)
    
    
    cd(outputsourcedir)
    saveID = [subj, '_IC:', num2str(nIC(i))];
    
    source1 = [];
    source1.pow = source.avg.pow(:,nIC(i));
    
    s_in = source1; mask = 'pow'; thre = 0.3;
    %     Run_surfacemap_template
    
    s_in.pos     = template_grid.pos;
    s_in.dim     = template_grid.dim;
    s_in.inside  = template_grid.inside;
    
    %%
    cfg = [];
    cfg.parameter = mask;
    % cfg.interpmethod = 'sphere_avg';
    cfg.interpmethod = 'smudge';
    cfg.coordsys     = 'mni';
    source_int = ft_sourceinterpolate(cfg, s_in, template_mri);
    
    cfg = [];
    %     cfg.subj = [];
    cfg.subj = saveID;
    cfg.mask = mask;
    cfg.thre = thre;
    if sIC==1
        cfg.saveflag = 1;
    else
        cfg.saveflag = 2;
    end
    cfg.colorbar = 2;
    vy_mapvisualisation(cfg, source_int);
    
    
end

%% mean value
saveID = [subj, '_mean_IC'];

source1 = [];
source1.pow = mean(source.avg.pow(:,nIC),2);

s_in = source1; mask = 'pow'; thre = 0.3;
%     Run_surfacemap_template

s_in.pos     = template_grid.pos;
s_in.dim     = template_grid.dim;
s_in.inside  = template_grid.inside;

%
cfg = [];
cfg.parameter = mask;
% cfg.interpmethod = 'sphere_avg';
cfg.interpmethod = 'smudge';
cfg.coordsys     = 'mni';
source_int = ft_sourceinterpolate(cfg, s_in, template_mri);

cfg = [];
%     cfg.subj = [];
cfg.subj = saveID;
cfg.mask = mask;
cfg.thre = thre;
if sIC==1
    cfg.saveflag = 1;
else
    cfg.saveflag = 2;
end
cfg.colorbar = 2;
vy_mapvisualisation(cfg, source_int);

%%
pause(2)

