function vy_source_dics_stats(cfg_main, ep_data)

cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.foilim = [2 40];
cfg.plotflag  = 2;
cfg.tapsmofrq = 1;
cfg.taper    = 'hanning';
cfg.output = 'fourier';
f_data.bsl = vy_fft(cfg, ep_data.bsl); f_data.bsl.elec = cfg_main.sens;
f_data.pst = vy_fft(cfg, ep_data.pst); f_data.pst.elec = cfg_main.sens;

% % freq analysis - prepration for DICS source analysis
% f_data.bsl = vy_fft(ep_data.bsl, [2,40], 0,[],0); f_data.bsl.elec = sens;
% f_data.pst = vy_fft(ep_data.pst, [2,40], 0,[],0); f_data.pst.elec = sens;

% PSD - sensor space
psd_bsl = squeeze(mean(mean(abs(f_data.bsl.fourierspctrm),2),1));
psd_pst = squeeze(mean(mean(abs(f_data.pst.fourierspctrm),2),1));
ff = linspace(1, 40, length(psd_pst));

figure,plot(ff,psd_bsl)
hold on
% plot(ff,psd,'g')
plot(ff,psd_pst,'r')
hold on
plot(ff,psd_pst - psd_bsl,'k')
xlabel('Hz'); ylabel('psd'),legend({'bsl','pst','diff'})

% outputdir_dics = fullfile(cfg_main.outputdir,'dics');
outputdir_dics = cfg_main.outputdir;
if exist(outputdir_dics, 'file') == 0, mkdir(outputdir_dics), end

% savepath = fullfile(outputdir_dics,['psd_',subj,'_',run,'.mat']);
% save(savepath, 'ff','psd_bsl','psd_pst', '-v7.3');
% savepath = fullfile(outputdir_dics,['fftcon_',subj,'_',run]);
% hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);

[a,b] = min(psd_pst - psd_bsl);
% f = ff(b); f = round(f);

% f = input('FOI? ');
f = 22;
cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.foilim = [f f];
cfg.plotflag  = 2;
cfg.tapsmofrq = 4;
cfg.taper    = 'dpss';
cfg.output = 'fourier';
% cfg.taper    = 'hanning';
[f_data.app,~,tapsmofrq] = vy_fft(cfg, ep_data.app); f_data.app.elec = cfg_main.sens;
f_data.bsl = vy_fft(cfg, ep_data.bsl); f_data.bsl.elec = cfg_main.sens;
f_data.pst = vy_fft(cfg, ep_data.pst); f_data.pst.elec = cfg_main.sens;

% f = 20;
% % Freq of interest - prepration for DICS source analysis
% f_data.app = vy_fft(ep_data.app, [f,f], 0,[],0); f_data.app.elec = sens;
% f_data.bsl = vy_fft(ep_data.bsl, [f,f], 0,[],0); f_data.bsl.elec = sens;
% f_data.pst = vy_fft(ep_data.pst, [f,f], 0,[],0); f_data.pst.elec = sens;

%%
cfg = [];
cfg.headmodel = cfg_main.headmodel;
% cfg.sourcemodel = cfg_main.sourcemodel;
cfg.grid = cfg_main.grid;
cfg.mtag = cfg_main.mtag;
s_data_dics = vy_source_freq(cfg, f_data);
% [s_data_dics, ~] = vy_source_freq(f_data, cfg_main.grid, cfg_main.headmodel, 'dics_stat');
% stat = vy_source_stat_montcarlo(s_data_dics);


s_data_dics.bsl.bnd = cfg_main.headmodel.bnd;
s_data_dics.pst.bnd = cfg_main.headmodel.bnd;

s_data_dics.pst.dim = cfg_main.grid.dim;
s_data_dics.pst.dim = cfg_main.grid.dim;
stat = vy_source_stat_montcarlo_cluster(s_data_dics);

%%
% cfg = []
% cfg.method =  'triangulation';
% 
% neighbours = ft_prepare_neighbours(cfg, s_data_dics.pst);
% 
% %%
% s = [];
% s = s_data_dics.pst; % timePost is the output of ft_timelockanalysis.
% s.grad.chanpos = cfg_main.grid.pos(cfg_main.grid.inside,:);
% s.grad.chanori = s.grad.chanpos;
% % s.grad.chanunit(1:size(s.grad.chanpos,1)) = s.grad.chanunit(1);
% % s.grad.chantype(1:size(s.grad.chanpos,1)) = s.grad.chantype(1);
% s.grad.tra    = cfg_main.headmodel.bnd;
% for i=1:size(s.grad.chanpos,1)
%     s.grad.label{i} = num2str(i);
%     s.label{i} = num2str(i);
%     s.grad.labelold =  num2str(i);
% end
% % prepare_neighbours determines what sensors may form clusters
% cfg_neighb.method   = 'distance';
% % cfg.method      = 'triangulation';
% neighbours          = ft_prepare_neighbours(cfg_neighb, s);

%%

% ft_neighbourplot([], stat)

stat.pos     = cfg_main.template_grid.pos;
stat.dim     = cfg_main.template_grid.dim;
stat.inside  = cfg_main.template_grid.inside;

tmp = stat.stat;
tmp2 = zeros(size(stat.pos,1),1);
tmp2(stat.inside) = tmp;

stats1  = stat;
stats1.stat =  tmp2;
stats1.mask = stat.inside;

% disp('1: Postive effects')
% disp('2: Negative effects')
% disp('3: Both effects')
% effect = input(':');
% 
stats2 = stats1;
% if effect == 1
%     stats2.stat(stats1.stat<0)=0;
% elseif effect == 2
%     stats2.stat(stats1.stat>0)=0;
% elseif effect == 3
%     stats2.stat = stats1.stat;
% end


stats2.stat(stats2.stat>0)=0;
stats2.stat(isnan(stats2.stat))=0;

%%
% idx = [];
% idx = find(stats2.prob > 0.025);
% stats3 = stats2;
% stats3.stat(idx) = 0;

%%
% param = [];
% param.mask = 'stat';
% param.loc = 'min';
% source_diff_dics = vy_source_plot(stats1,cfg_main.template_mri,param,2);

%%
% cfg = [];
% cfg.parameter = 'pow';
% cfg.operation = '(x1-x2)/(x1+x2)';
% source_diff_dics = ft_math(cfg,s_data_dics.pst,s_data_dics.bsl);
% source_diff_dics.pos     = cfg_main.template_grid.pos;
% source_diff_dics.dim     = cfg_main.template_grid.dim;
% source_diff_dics.inside  = cfg_main.template_grid.inside;
% source_diff_dics.pow(source_diff_dics.pow>0)=0;


%% save group average as nii

% addpath(allpath.connpath);
% addpath(allpath.spm_path);
% close all
% 
% projthresh = 0.6;
% source_dics1 = vy_vol_thresh(source_int_dics,projthresh,'stat'); % abs
% source_dics1.stat = -(source_dics1.stat./max(source_dics1.stat(:)));
% 
% savenii = ['test.nii'];
% vy_savenifti(source_dics1,'stat',savenii);
% 
% Opt = [];
% Opt.savenii = 0; Opt.savefig = 0;
% Opt.savename = ['Test'];
% vy_surfce_vis2(source_int_dics,['test.nii'], Opt);
% vy_surfce_vis(source_int_dics,'test.nii', Opt)
% 
% 
% %-Surf-vis
% opt = [];
% opt.run = [];
% opt.tsk = [];
% opt.subj = '1';
% opt.savedir = [];
% opt.savenii = 2;
% %     opt.plot = '-mosaic';
% opt.plot = '-row';
% vy_surfce_vis([],'test.nii', opt);

%%
toi = cfg_main.toi;
if ~isempty(cfg_main.flag.savetag)
    save([cfg_main.savedata,'_',cfg_main.flag.savetag '_', num2str(f),'Hz.mat'], 'stats2', 'f','toi','-v7.3');
else
    save([cfg_main.savedata,'_',num2str(f),'Hz.mat'], 'stats2', 'f','toi','-v7.3');
end

%     savefig = fullfile(outputdir_dics,[num2str(f-tapsmofrq),'_',num2str(f+tapsmofrq),'Hz','_1_',cfg_main.subj]);
if size(toi,1) == 1
    savefig = fullfile(outputdir_dics,['stat_',num2str(toi(1)),'_',num2str(toi(2)),'sec_',num2str(f),'Hz','_1_',cfg_main.subj]);
else
    if ~isempty(cfg_main.flag.savetag)
        savefig = fullfile(outputdir_dics,['stat_', num2str(toi(2,1)),'_',num2str(toi(2,2)),'sec_',num2str(f),'Hz_',cfg_main.flag.savetag,'_1_',cfg_main.subj]);
    else
        savefig = fullfile(outputdir_dics,['stat_',num2str(toi(2,1)),'_',num2str(toi(2,2)),'sec_',num2str(f),'Hz','_1_',cfg_main.subj]);
    end
end
% outputdir_dics = fullfile(outputdir,'dics');
% if exist(outputdir_dics, 'file') == 0, mkdir(outputdir_dics), end
% savedata = fullfile(outputdir_dics,['s_dics_',subj,'_',run,'.mat']);
% save(outputdir_dics, 'source_diff_dics', '-v7.3');

% mtd = 'source_dics_stats';
% savefig = fullfile(outputdir_dics,[num2str(f),'Hz','_1_',cfg_main.subj]);
% savefig = fullfile(outputdir_dics,[num2str(f-tapsmofrq),'_',num2str(f+tapsmofrq),'Hz','_1_',cfg_main.subj]);
% savefig = fullfile(outputdir_dics,[num2str(cfg_main.toi(2,1)),'_',num2str(cfg_main.toi(2,2)),'sec_',num2str(f-tapsmofrq),'_',num2str(f+tapsmofrq),'Hz','_1_',cfg_main.subj]);


cfg = [];
cfg.mask = 'stat';
cfg.loc  = 'min';
cfg.template = cfg_main.template_mri;
cfg.savefile = savefig;
cfg.volnorm  = 2; % yes: 1
source_int_dics = vy_source_plot(cfg, stats2);
pause(1),

% param = [];
% param.mask = 'stat';
% param.loc = 'min';
% source_int_dics = vy_source_plot(stats1,cfg_main.template_mri,param,2);
% savefig = fullfile(outputdir_dics,[num2str(f),'Hz','_1_',cfg_main.subj]);
% % hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);
% print(savefig,'-depsc')

clear savepath
savepath{1} = fullfile(outputdir_dics,['stat_', num2str(cfg_main.toi(2,1)),'_',num2str(cfg_main.toi(2,2)),'sec_',num2str(f),'Hz','_2_',cfg_main.subj]);
savepath{2} = fullfile(outputdir_dics,['stat_', num2str(cfg_main.toi(2,1)),'_',num2str(cfg_main.toi(2,2)),'sec_',num2str(f),'Hz','_3_',cfg_main.subj]);

cfg = [];
cfg.subj = cfg_main.subj;
cfg.mask = 'stat';
cfg.thre = 0.6;
cfg.savepath = savepath;
vy_mapvisualisation(cfg, source_int_dics);
% vy_mapvisualisation_light(cfg, source_int_dics);

% vy_mapvisualisation(source_int_dics,cfg.mask,0.6, savepath);
% vy_mapvisualisation(source_int_dics,cfg.mask,0.6, []);

% restoredefaultpath
% addpath((cfg_main.allpath.ft_path));
% ft_defaults
% addpath(genpath(cfg_main.allpath.hcp_path));
% addpath(genpath(cfg_main.allpath.cd_org));
% addpath(genpath(cfg_main.allpath.exfig_path));
%
% cfg = [];
% cfg.maskparam = 'pow';
% cfg.save.savepath =  savepath;
% % cfg.saveformat = '-eps';
% cfg.save.saveformat = '-png';
% cfg.save.pixdim     = 12;
% cfg.projthresh      = 0.6;
% vy_surfmap(cfg, source_int_dics);




