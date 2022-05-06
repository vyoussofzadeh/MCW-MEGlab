
outputdir1 = fullfile(outd.sub, 'spm_source');
if exist(outputdir1, 'file') == 0
    mkdir(outputdir1);   % create a directory
end
addpath(genpath(allpath.spm_path))
cfg = [];
cfg.toilim = [-0.4 1.2];
eint_data = ft_redefinetrial(cfg, cln_data);

%%
if exist(mripfile, 'file') == 2
    cd(outputdir1);
    cfg = [];
    cfg.allpath = allpath;
    cfg.datafile = datafile;
    cfg.eint_data = eint_data;
    cfg.sub       = subj;
    cfg.mripfile = mripfile;
    vy_forward_spm_meg(cfg);
    cd(outputdir1)
end
cd ..

% revert to the newer ft!
restoredefaultpath
addpath((allpath.ft_path));
ft_defaults
addpath(genpath(allpath.hcp_path));
addpath(genpath(allpath.cd_org));

%% headmodel from SPM (surface-based)
spm_path = allpath.spm_path;
disp('Source model spm2ft');
vy_spm2ft_headmodel

%%
figure; hold on;
ft_plot_vol(individual_headmodel, 'facecolor', 'none'); alpha 0.5;
ft_plot_mesh(sourcemodel, 'facecolor', 'cortex', 'edgecolor', 'none'); camlight;
hold on;
ft_plot_headshape(headshape);

%%
outputdir1 = fullfile(outd.sub, 'Surface_source');
if exist(outputdir1, 'file') == 0, mkdir(outputdir1), end
mtd = 'dics';
d = ['.\spm_source\m',subj];
D = spm_eeg_load(d);
allmeshvert_mni = D.inv{1}.mesh.tess_mni.vert;
cd(outputdir1);


%% Spectral analysis
datain = cln_data;
cfg_main = [];
cfg_main.fmax = fmax;
cfg_main.sens = datain.grad;
cfg_main.outputdir = savepath;
cfg_main.freq_of_interest  = freq_of_interest; % Hz

cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.foilim = [2 cfg_main.fmax];
cfg.plotflag  = 2;
cfg.tapsmofrq       = 1;
cfg.taper    = 'hanning';
f_data.bsl = vy_fft(cfg, ep_data.bsl); f_data.bsl.elec = cfg_main.sens;
f_data.pst = vy_fft(cfg, ep_data.pst); f_data.pst.elec = cfg_main.sens;

% PSD - sensor space
psd_bsl = squeeze(mean(mean(abs(f_data.bsl.fourierspctrm),2),1));
psd_pst = squeeze(mean(mean(abs(f_data.pst.fourierspctrm),2),1));
ff = linspace(1, cfg_main.fmax, length(psd_pst));


outputdir_dics = cfg_main.outputdir;
% if exist(outputdir_dics, 'file') == 0, mkdir(outputdir_dics), end

f_sugg = round(cfg_main.freq_of_interest);
disp(['Suggested by TFR: ', num2str(f_sugg),'(+-3Hz)']);
disp(['Select foi,eg ,',num2str(f_sugg),':']);
clear('input')
f = input('Freq of interest? ');
%         tapsmofrq = 4;
tapsmofrq = input('tapsmofrq, e.g. 4 Hz? ');

cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.foilim = [f f];
cfg.plotflag  = 2;
cfg.taper    = 'dpss'; cfg.tapsmofrq  = tapsmofrq;

if f < 4, cfg.tapsmofrq  = 1; cfg.taper    = 'hanning'; end

[f_data.app,~,~,~] = vy_fft(cfg, ep_data.app); f_data.app.elec = cfg_main.sens;
f_data.bsl = vy_fft(cfg, ep_data.bsl); f_data.bsl.elec = cfg_main.sens;
f_data.pst = vy_fft(cfg, ep_data.pst); f_data.pst.elec = cfg_main.sens;

%%
%- FFT_based
cfg = [];
cfg.method = 'dics';
cfg.dics.lambda = '100%';
cfg.sourcemodel  = individual_grid;
cfg.frequency    = f_data.app.freq;
cfg.headmodel = individual_headmodel;
cfg.dics.keepfilter = 'yes';
cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
sourceavg = ft_sourceanalysis(cfg, f_data.app);

cfg = [];
cfg.method = 'dics';
cfg.dics.lambda = '0%';
cfg.sourcemodel        = individual_grid;
cfg.sourcemodel.filter = sourceavg.avg.filter;
cfg.dics.fixedori    = 'yes'; % project on axis of most variance using SVD
cfg.headmodel = individual_headmodel;
s_data.bsl      = ft_sourceanalysis(cfg, f_data.bsl);
s_data.pst      = ft_sourceanalysis(cfg, f_data.pst);
%                 disp('desynchronisation effects: 1, synchronisation effects: 2:');
%                 ask_sd = input(':');
ask_sd = 1;


cfg = [];
cfg.parameter = 'pow';
cfg.operation = 'log10(x1/x2)'; % sourceA divided by sourceB
source_diff_dics = ft_math(cfg,s_data.pst,s_data.bsl);
% source_diff_dics.pow(isnan(source_diff_dics.pow))=0;
% source_diff_dics.pow(source_diff_dics.pow>0)=0;
% source_diff_dics.pow = abs(source_diff_dics.pow);

figure
m = source_diff_dics.pow;
bnd.pnt = sourcemodel.pos;
bnd.tri = sourcemodel.tri;
ft_plot_mesh(bnd, 'vertexcolor', abs(m));
colorbar


Jabs_DICS{val} = abs(source_diff_dics.pow);

[~,ii] = sort(Jabs_DICS{val},'descend');
Ndip = 250;
ii = ii(1:Ndip);
Is = D.inv{1}.inverse.Is;
ivert_dics{val} = Is(ii);

 Jabs = Jabs_DICS{val}; in = ivert_dics{val};
spm_mip(Jabs(in),allmeshvert_mni(in,:)',6); axis image; colormap gray;
   


