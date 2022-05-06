function cln_data = SM_preprocess(data,lay,subj, allpath)


cfg = [];
cfg.hpfilter = 'yes';
cfg.lpfilter = 'yes';
cfg.dftfilter = 'yes';
cfg.hpfiltord = 3;
cfg.hpfreq = 1;
cfg.lpfreq = 40;
f_data = ft_preprocessing(cfg, data);

%%
% cfg = [];
% cfg.pflag    = 0; % yes:1, No:2
% cfg.saveflag = 0; % yes:1, No:2
% cfg.savepath = [];
% cfg.latency  = [f_data.time{1}(1),f_data.time{1}(end)];%[-200,900]./1000;
% cfg.rejectpercentage = .95;
% [r_data,~] = vy_artifactreject(cfg, f_data);
% % save('ica/ICAreport.mat','report_ic');
% % disp(report_ic)
r_data = f_data;

%%
cfg = [];
cfg.lay = lay;
cfg.subj = subj;
cfg.n = 20;
cfg.savefig = 0;
cfg.allpath = allpath;
comp = vy_ica(cfg, r_data);

bic = input(['Select bad ICs for ' subj,':']);
cfg.component = comp.label(bic);
cfg.updatesens = 'no';
% cfg.trials = trials;
cln_data = ft_rejectcomponent(cfg, comp, r_data);
close all

%% back to new ft
restoredefaultpath
addpath((allpath.ft_path));
ft_defaults
addpath(genpath(allpath.hcp_path));
addpath(genpath(allpath.cd_org));
addpath(genpath(allpath.exfig_path));

