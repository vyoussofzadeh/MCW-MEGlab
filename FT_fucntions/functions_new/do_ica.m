function dataout = do_ica(cfg_main, data)

cfg            = [];
cfg.method     = 'runica';
cfg.numcomponent = cfg_main.n;       % specify the component(s) that should be plotted
comp           = ft_componentanalysis(cfg, data);

%%
restoredefaultpath
addpath(genpath(cfg_main.allpath.ft_path18));
% addpath(genpath(cfg_main.allpath.hcp_path));
addpath(genpath(cfg_main.allpath.script_path));

cfg = [];
cfg.viewmode = 'component';
cfg.layout = cfg_main.lay;
ft_databrowser(cfg, comp);
colormap(brewermap(256, '*RdYlBu'));
set(gcf, 'Position', [600   600   700   500]);
set(gcf,'name',cfg_main.subj,'numbertitle','off')

%%
cfg              = [];
cfg.output       = 'pow';
cfg.channel      = 'all';%compute the power spectrum in all ICs
cfg.method       = 'mtmfft';
cfg.taper        = 'hanning';
cfg.foi          = 2:2:30;
freq = ft_freqanalysis(cfg, comp);

%%
n = cfg_main.n;
nby1 = 5; nby2 = 4;

Nfigs = ceil(size(comp.topo,1)/n);
tot = Nfigs*n;

rptvect = 1:size(comp.topo,1);
rptvect = padarray(rptvect, [0 tot-size(comp.topo,1)], 0,'post');
rptvect = reshape(rptvect,n,Nfigs)';

figure
for r=1:n
    cfg=[];
    cfg.channel = rptvect(:,r);
    subplot(nby1,nby2,r); set(gca,'color','none');
    ft_singleplotER(cfg,freq);
end
colormap(brewermap(256, '*RdYlBu'));
set(gcf, 'Position', [800   600   800   500]);
set(gcf,'name',cfg_main.subj,'numbertitle','off')


% Remove bad ICA components
cfg = [];
cfg.updatesens = 'no';
bic = input(['Select bad ICs for ',num2str(cfg_main.subj),':']);
cfg.component = comp.label(bic);
dataout = ft_rejectcomponent(cfg, comp, data);

end
