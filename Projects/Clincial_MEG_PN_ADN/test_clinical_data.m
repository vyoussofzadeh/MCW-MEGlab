
clear,
restoredefaultpath
ft_path = '/MEG_data/MEG_Tools/fieldtrip/fieldtrip_2022';
addpath(ft_path)
ft_defaults


datafile = '/MEG_data/epilepsy/schubert_daniel/260504/sss/Run10_DFNM_supine/Run10_DFNM_supine_raw_t_sss.fif';
subj = 'schubert_daniel';

fg = []; cfg.layout = 'neuromag306mag.lay'; lay = ft_prepare_layout(cfg);
% ft_layoutplot(cfg);

addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/FT_fucntions/functions');
addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/FT_fucntions/run');
addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/Projects/psb_pilot/functions/ft_private')
addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/FT_fucntions/helper')

% outd.sub = '/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/Projects/Clincial_MEG_PN_ADN';
% subj = 'sd';
% ic_selection = 1; % 1: manual, 2: automated
% flag.preprocessing.filtering = 1;
% Run_preprocess

cd('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/Projects/Clincial_MEG_PN_ADN/')

% cfg = [];
% cfg.eventid = 5;
% cfg.epochtype = 'STI001';
% cfg.datafile  = datafile;
% cfg.hpfreq = 0.1;
% cfg.lpfreq = 40;
% [f_data] = vy_preprocess(cfg);
%
% cfg = [];
% cfg.pflag = 1; % yes:1, No:2
% cfg.saveflag = 1; % yes:1, No:2
% cfg.savepath = savepath;
% cfg.latency = [f_data.time{1}(1),f_data.time{1}(end)];%[-200,900]./1000;
% cfg.rejectpercentage = .95;
% [r_data,report] = vy_artifactreject(cfg, f_data);
%
% cfg = [];
% cfg.savepath = outd.sub;
% cfg.savefile = fullfile(outd.sub,['ic_',subj,'.mat']);
% cfg.saveflag = 1;
% cfg.overwrite = 2;
% cfg.lay = lay;
% cfg.n   = 20;
% cfg.subj = subj;
% cfg.allpath = [];
% cfg.select = ic_selection;
% ic_data = vy_ica_cleaning_light2(cfg, r_data);

load('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/Projects/Clincial_MEG_PN_ADN/ic_sd.mat');

dtag = []; dtag.dcon = 'dfnm';
fmax = 40;
datain = cln_data;
% flag.freq = 1;
% if flag.freq == 1
%     stag = 'tsk_baseline';
%     Run_freq
%     disp(['time_of_interest:',num2str(time_of_interest),'sec']);
%     disp(['freq_of_interest:',num2str(freq_of_interest),'Hz']);
%     L = 0.3;
% end
%
% task = 1;
%
% % -Time-locked
% flag.time = 1;
% if flag.time == 1
%     %--setting baseline interval
%     toi(1,:) = [-0.3,0];
%
%     %-- setting the post-stim interval
%     disp(['suggested time interval:',num2str(time_of_interest), '+-', num2str(L),' Sec']);
%     %     disp('Yes: 1, No: 2');
%     %     tfa = input('Is it OK to proceed?');
%     %     disp('the following time was selected');
%     %        if (time_of_interest-L) > 0.4 && (time_of_interest+L) < 2.3
%     if (time_of_interest-L) >= 0.3 && (time_of_interest+L) < 2.3
%         toi(2,:) = round([time_of_interest-L, time_of_interest+L],1);
%     else
%         switch task
%             case 1
%                 %                         toi = [-0.3,0;1.1,1.7]; % Best of DN
%                 toi(2,:) = [0.8,1.5]; % Best of DN
%             case 2
%                 toi(2,:) = [0.4,1.2]; % Best for PN, left IFG
%                 toi(2,:) = [0.4,1.6]; % Best for PN, left IFG
%                 toi(2,:) = [0.7,1.6]; % Best for PN, left IFG
%                 toi(2,:) = [1,1.6]; % Best for PN, left IFG
%         end
%     end
%     disp(['[',num2str(toi(1,1)),',',num2str(toi(1,2)),'] sec interval was selected as bsl']);
%     disp(['[',num2str(toi(2,1)),',',num2str(toi(2,2)),'] sec interval was selected as pst']);
%     Run_time
%
% end
%
% Run_grandmean

%% Anatomy
restoredefaultpath
ft_path = '/MEG_data/MEG_Tools/fieldtrip/fieldtrip_2022';
addpath(ft_path)
ft_defaults
addpath(genpath('/MEG_data/MEG_Tools/fieldtrip/fieldtrip_2022/plotting'))
addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git_BAK/Tools/BS/external/')
addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/Projects/CIDMEG_BGross/FT pipelines/run');
addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/Projects/CIDMEG_BGross/FT pipelines/functions');
addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/Projects/Clincial_MEG_PN_ADN/func')
atlas = '/MEG_data/MEG_Tools/fieldtrip/fieldtrip_2022/template/atlas';
mridir = '/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/Projects/Clincial_MEG_PN_ADN/nii';
mrifile = 'mri_resliced.nii';
data_clean = datain;
outdir = '/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/Projects/Clincial_MEG_PN_ADN/Dan_sh';
addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/FT/Data_file')
Run_anatomy

%%
addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/FT_fucntions/functions');
outputmridir = '/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/Projects/Clincial_MEG_PN_ADN/Dan_sh/schubert_daniel/anat';
cfg = [];
cfg.saveflag = 0;
cfg.headmodel = individual_headmodel;
cfg.leadfield = individual_grid;
cfg.mri_realigned  = mri_realigned;
cfg.headshape = headshape;
cfg.outputmridir = outputmridir;
cfg.mtd = 'vol';
vy_mri_inspection(cfg, data_clean);

%%
template_mri = ft_read_mri(fullfile(ft_path,'template/anatomy','single_subj_T1.nii')); %

%%
sens = ft_read_sens(datafile);
sens = ft_convert_units(sens,'mm');

%%
toi = [-0.3,0;0.3,0.9];
% toi = [-0.3,0;0.3, 0.8];
ep_data = vy_epoch(data_clean, toi);

cfg = [];
ep_data.app = ft_appenddata(cfg,ep_data.bsl,ep_data.pst);

%%
freq_of_interest = 20;
flag.savetag = 0;
mtag = 'dics';
toi = [0.3,0.8];
cfg = [];
cfg.grid = individual_grid;
cfg.allpath = [];
cfg.freq_of_interest  = freq_of_interest; % Hz
cfg.headmodel = individual_headmodel;
cfg.sens = sens;
cfg.mtag = mtag;
cfg.subj = subj;
cfg.toi = toi;
cfg.outputdir = '/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/Projects/Clincial_MEG_PN_ADN/Dan_sh';
switch choose_grid
    case 1
        cfg.template_grid = [];
        cfg.mri_aligned = outanat.mri_realigned;
    case 2
        cfg.template_grid = template_grid;
end
cfg.template_mri = template_mri;
cfg.fmax = fmax;
cfg.savedata = fullfile(cfg.outputdir,[mtag,'_',subj]);
cfg.flag = flag;
cfg.plotting =1;
vy_source_dics(cfg, ep_data);

%%
% meshgridres = 1;
% 
% indir = '/MEG_data/epilepsy/';
% subj = 'schubert_daniel';
% outputmridir = '/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/Projects/Clincial_MEG_PN_ADN';
% 
% fmax = 40;
% % Volumetric-based analysis
% anatomy_check_flag = 1;
% Run_volumetric

