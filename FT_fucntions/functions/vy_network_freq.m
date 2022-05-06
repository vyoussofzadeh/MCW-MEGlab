function vy_network_freq(cfg_main, ep_data)

Run_fft_4dics

%%
% cfg = [];
% cfg.headmodel = cfg_main.headmodel;
% % cfg.sourcemodel = cfg_main.sourcemodel;
% cfg.grid = cfg_main.grid;
% cfg.mtag = cfg_main.mtag;
% cfg.kappa = [];
% s_data_dics = vy_source_freq(cfg, f_data);
% % s_data_dics = vy_source_freq(f_data, cfg_main.grid, cfg_main.headmodel, cfg_main.mtag);

mtd = 'pcc';
%- FFT_based
cfg = [];
cfg.method = mtd;
cfg.(mtd).lambda = '10%';
cfg.senstype  = 'MEG';
cfg.sourcemodel  = cfg_main.grid;
cfg.frequency = f_data.pst.freq;
cfg.headmodel = cfg_main.headmodel;
%                 cfg.dics.keepfilter = 'yes';
cfg.(mtd).fixedori = 'yes'; % project on axis of most variance using SVD
%                 cfg.dics.projectnoise = 'yes';
s_data.bsl      = ft_sourceanalysis(cfg, f_data.bsl);
s_data.pst      = ft_sourceanalysis(cfg, f_data.pst);

% cfg = [];
% cfg.method = mtd;
% cfg.pcc.lambda = '10%';
% cfg.frequency    = f_data.app.freq;
% cfg.sourcemodel = cfg_main.grid;
% cfg.headmodel = cfg_main.headmodel;
% cfg.(mtd).keepfilter = 'yes';
% cfg.(mtd).fixedori    = 'yes'; % project on axis of most variance using SVD
% %         cfg.kappa = cfg_main.kappa;
% sourceavg = ft_sourceanalysis(cfg, f_data.app);
% 
% cfg = [];
% cfg.method = 'pcc';
% %         cfg.dics.lambda = '5%';
% cfg.sourcemodel = cfg_main.grid;
% cfg.sourcemodel.filter = sourceavg.avg.filter;
% % cfg.(mtd).fixedori    = 'yes'; % project on axis of most variance using SVD
% cfg.headmodel = cfg_main.headmodel;
% %         cfg.kappa = cfg_main.kappa;
% s_data.bsl      = ft_sourceanalysis(cfg, f_data.bsl);
% s_data.pst      = ft_sourceanalysis(cfg, f_data.pst);

%%
disp('1:coh');
disp('2:imgcoh');
disp('3:plv');
disp('4:powcorr');
disp('5:powcorr_ortho');
disp('6:wpli_debiased');
disp('7:wpli');
conn_sel = input(':');
cfg           = [];
switch conn_sel
    case 1
        cfg.method    = 'coh'; cfg.complex = 'abs'; label = 'coh'; savetag = 'coh';
    case 2
        cfg.method    = 'coh'; cfg.complex = 'absimag'; label = 'coh'; savetag = 'imcoh';
    case 3
        cfg.method    = 'plv'; cfg.complex = 'abs'; label = 'plv'; savetag = 'plv';
    case 4
        cfg.method    = 'powcorr'; cfg.complex = 'abs'; label = 'powcorr'; savetag = 'powcorr';
    case 5
        cfg.method    = 'powcorr_ortho'; cfg.complex = 'abs'; label = 'powcorr'; savetag = 'powcorr-ortho';
    case 6
        cfg.method    = 'wpli_debiased'; cfg.complex = 'abs'; label = 'wpli_debiased'; savetag = 'wpli-debiased';
    case 7
        cfg.method    = 'wpli'; cfg.complex = 'abs'; label = 'wpli'; savetag = 'wpli';
end
% cfg.method    = 'wpli'; cfg.complex = 'abs'; label = 'wpli';
% cfg.method    = 'psi'; cfg.complex = 'abs'; label = 'psi'; cfg.bandwidth = 4;
% cfg.method    = 'wppc'; cfg.complex = 'abs'; label = 'wppc';


tmp = cat(1, s_data.bsl.avg.mom{s_data.bsl.inside});
f_data1.bsl = f_data.bsl;
f_data1.bsl.fourierspctrm = tmp';
lbl = [];
for i=1:size(f_data1.bsl.fourierspctrm,2)
    lbl{i}=num2str(i);
end
f_data1.bsl.label = lbl';


tmp = cat(1, s_data.pst.avg.mom{s_data.pst.inside});
f_data1.pst = f_data.pst;
f_data1.pst.fourierspctrm = tmp';
lbl = [];
for i=1:size(f_data1.pst.fourierspctrm,2)
    lbl{i}=num2str(i);
end
f_data1.pst.label = lbl';

conn_bsl           = ft_connectivityanalysis(cfg, f_data1.bsl); conn_bsl.dimord    = 'pos_pos';
conn_pst           = ft_connectivityanalysis(cfg, f_data1.pst); conn_pst.dimord    = 'pos_pos';

%                 conn_bsl           = ft_connectivityanalysis(cfg, s_data.bsl); conn_bsl.dimord    = 'pos_pos';
%                 conn_pst           = ft_connectivityanalysis(cfg, s_data.pst); conn_pst.dimord    = 'pos_pos';

%%
mtd = [label,'spctrm']; % method par.
% gtm = 'eigenvector_cent'; % grapth theory measure
gtm = 'degrees';

figure,imagesc(conn_pst.(mtd) - conn_bsl.(mtd)); colorbar
%                 figure,imagesc(conn_bsl.(mtd)); colorbar
%                 figure,imagesc(conn_pst.(mtd)); colorbar

%                 conn_bsl.(mtd)(isnan(conn_bsl.(mtd))) = 0;
%                 conn_pst.(mtd)(isnan(conn_pst.(mtd))) = 0;
%
%                 tmp = conn_bsl.(mtd); tmp(isnan(tmp(:)))=0; tmp = tmp./max(tmp(:)); conn_bsl.(mtd) = tmp;
%                 tmp = s_data.pst.avg.pow; tmp(isnan(tmp))=0; tmp = tmp./max(tmp(:)); s_data.pst.avg.pow = tmp;


tmp = conn_bsl.(mtd); tmp = tmp./max(tmp(:)); tmp(isnan(tmp(:))) = 0;conn_bsl.(mtd) = tmp;
tmp = conn_pst.(mtd); tmp = tmp./max(tmp(:)); tmp(isnan(tmp(:))) = 0; conn_pst.(mtd) = tmp;

thre = 0.5;
cfg = [];
cfg.method    = gtm;
cfg.parameter = mtd;
cfg.threshold = thre.*max(conn_bsl.(mtd)(:));
net_bsl = ft_networkanalysis(cfg,conn_bsl);
net_bsl.pos = s_data.bsl.pos;
net_bsl.dim  = s_data.bsl.dim;
net_bsl.inside  = s_data.bsl.inside;

cfg = [];
cfg.method    = gtm;
cfg.parameter = mtd;
cfg.threshold = thre.*max(conn_pst.(mtd)(:));
net_pst = ft_networkanalysis(cfg,conn_pst);
net_pst.pos = s_data.pst.pos;
net_pst.dim  = s_data.pst.dim;
net_pst.inside  = s_data.pst.inside;

%%
cfg = [];
cfg.parameter = gtm;
cfg.operation = 'x1-x2';
net_diff = ft_math(cfg,net_pst,net_bsl);

net_diff.(gtm) = zeros(size(net_diff.inside,1),1);
net_diff.(gtm)(net_diff.inside) = net_pst.(gtm) - net_bsl.(gtm);

%                 net_diff.(gtm) = net_pst.(gtm) - net_bsl.(gtm);
% net_diff.(gtm)(net_diff.(gtm)<0)=0;
% net_diff.(gtm)(net_diff.(gtm)>0)=0;
% net_diff = net_pst;
%                 net_diff.elec = freq_sel1.elec;

net_diff.pos     = cfg_main.template_grid.pos;
net_diff.dim     = cfg_main.template_grid.dim;
net_diff.inside  = cfg_main.template_grid.inside;
%                 net_diff.pow = (s_data.pst.avg.pow - s_data.bsl.avg.pow)./(s_data.pst.avg.pow + s_data.bsl.avg.pow);
net_diff.(gtm)(net_diff.(gtm)<0)=0;

%                 source_diff_dics.pow(source_diff_dics.pow>0)=0;
savefig = fullfile(cfg_main.outputdir,savetag);

cfg = [];
cfg.mask = gtm;
cfg.loc = 'min';
cfg.template = cfg_main.template_mri;
cfg.savefile = savefig;
cfg.volnorm     = 2; % yes: 1
net_diff1 = vy_source_plot(cfg, net_diff);

%%
% savepath = [];
% savepath{1} = 'conn2.png';
% savepath{2} = 'conn3.png';

clear savepath
savepath{1} = fullfile(cfg_main.outputdir,['R_',savetag]);
savepath{2} = fullfile(cfg_main.outputdir,['L_',savetag]);

cfg = [];
cfg.subj = cfg_main.subj;
cfg.mask = gtm;
cfg.thre = 0.6;
cfg.savepath = savepath;
cfg.colorbar = 2;
vy_mapvisualisation(cfg, net_diff1);

%%
% using older ver of ft for network analysis
% restoredefaultpath
% addpath(genpath(cfg_main.allpath.ft_oldp.ft_old));
% addpath(genpath(cfg_main.p.hcp_path));
% addpath(genpath(cfg_main.p.cd_org));
%
% restoredefaultpath
% addpath(genpath(cfg_main.allpath.ft_old));
% addpath(genpath(cfg_main.allpath.hcp_path));
% addpath(genpath(cfg_main.allpath.cd_org));

%%
% mtd = 'plv';
% mtd_par = 'plvspctrm';
% gtm = 'eigenvector_cent';
% 
% conn_par = [];
% conn_par.method   = mtd;
% conn_par.idx      = mtd_par;
% conn_par.complex  = [];
% 
% net_par.gtm       = gtm ; in2 = 2;
% net_par.threshold = 0.7;
% [source_conn_bsl, network_bsl] = vy_conn(s_data.bsl,conn_par,net_par);
% [source_conn_pst, network_pst] = vy_conn(s_data.pst,conn_par,net_par);
% 
% %%
% cfg = [];
% cfg.parameter = gtm;
% cfg.operation = 'x1-x2';
% network_diff_lcmv = ft_math(cfg,network_pst,network_bsl);
% network_diff_lcmv.pos     = cfg_main.template_grid.pos;
% network_diff_lcmv.dim     = cfg_main.template_grid.dim;
% network_diff_lcmv.inside  = cfg_main.template_grid.inside;
% network_diff_lcmv.(gtm)   = zscore(network_diff_lcmv.(gtm));
% network_diff_lcmv.eigenvector_cent(network_diff_lcmv.eigenvector_cent<0)=0;
% 
% %% Revert to new ft!
% restoredefaultpath
% addpath((cfg_main.allpath.ft_path));
% ft_defaults
% addpath(genpath(cfg_main.allpath.hcp_path));
% addpath(genpath(cfg_main.allpath.cd_org));
% addpath(genpath(cfg_main.allpath.exfig_path));
% 
% %%
% toi = cfg_main.toi;
% if size(toi,1) > 1
%     sname2 = [num2str(toi(2,1)),'_',num2str(toi(2,2)),'sec'];
%     disp([num2str(toi(2,1)),'_',num2str(toi(2,2)),'sec is analysing'])
% else
%     sname2 = [num2str(toi(1)),'_',num2str(toi(2)),'sec'];
%     disp([num2str(toi(1)),'_',num2str(toi(2)),'sec is analysing'])
% end
% 
% %%
% 
% savepath = fullfile(cfg_main.outputdir);
% if exist(savepath, 'file') == 0, mkdir(savepath), end
% savedata = fullfile(savepath,['n_',sname2,'_',cfg_main.subj,'.mat']);
% save(savedata, 'network_diff_lcmv', '-v7.3');
% 
% mtd = 'network_evc';
% 
% % oldmask = network_diff_lcmv.(gtm);
% % network_diff_lcmv.mask = (oldmask - min(oldmask(:))) ./ (max(oldmask(:)) - min(oldmask(:))); %
% % % set the new mask to range between 0 and 1
% % network_diff_lcmv.mask(oldmask < 0.9) = 0; % mask out non-significant voxels
% % savefig = fullfile(outputdir_dics,[num2str(f-tapsmofrq),'_',num2str(f+tapsmofrq),'Hz','_1_',cfg_main.subj]);
% % savefig = fullfile(savepath,[mtd,'_1_',cfg_main.subj]);
% savefig = fullfile(savepath,['n_',sname2,'_',cfg_main.subj]);
% 
% cfg = [];
% cfg.mask = gtm;
% cfg.loc = 'max';
% cfg.template = cfg_main.template_mri;
% cfg.savefile = savefig;
% cfg.volnorm     = 2; % yes: 1
% network_int = vy_source_plot(cfg, network_diff_lcmv);
% 
% % vy_mapvisualisation(network_int_lcmv,gtm,0.6,savep)
% % vy_mapvisualisation(network_int,gtm,0.6, []);
% 
% clear savepath
% savepath{1} = fullfile(cfg_main.outputdir,['_R_',cfg_main.subj]);
% savepath{2} = fullfile(cfg_main.outputdir,['_L_',cfg_main.subj]);
% % vy_mapvisualisation(network_int,gtm,0.6, savepath);
% % vy_mapvisualisation(network_int,gtm, 0.6, []);
% 
% cfg = [];
% cfg.subj = cfg_main.subj;
% cfg.mask = gtm;
% cfg.thre = 0.6;
% cfg.savepath = savepath;
% cfg.savepath = [];
% vy_mapvisualisation(cfg, network_int);
% 
% % savenii = fullfile(savepath,['n_',subj,'.nii']);
% % vy_savenifti(network_int_lcmv, gtm, savenii);

%%


% param = [];
% param.mask = gtm;
% param.loc = 'max';
% network_int_lcmv = vy_source_plot(network_diff_lcmv,cfg_main.template_mri,param,2);
% hcp_write_figure([savefig,'.png'], gcf, 'resolution', 300);
% saveformat = '-png';
% pixdim     = '-m8';
% export_fig(savefig, saveformat, pixdim)
% print(gcf, '-dpdf',[savefig,'.pdf']);

% clear savep
% savep{1} = fullfile(savepath,[mtd,'_2_',cfg_main.subj]);
% savep{2} = fullfile(savepath,[mtd,'_3_',cfg_main.subj]);

% cfg = [];
% cfg.maskparam = gtm;
% cfg.save.savepath =  savep;
% % cfg.saveformat = '-eps';
% cfg.save.saveformat = '-png';
% cfg.save.pixdim     = 12;
% cfg.projthresh      = 0.6;
% vy_surfmap(cfg, network_int_lcmv);
% cfg = [];
% cfg.funparameter = param.mask;
% cfg.method = 'ortho';  % plot slices
% ft_sourceplot(cfg, network_diff_lcmv);


%%
% network_int1 = network_int_lcmv;
% thre = 0.5;
% network_int1.eigenvector_cent(network_int1.eigenvector_cent < thre*max(network_int1.eigenvector_cent(:)))=0;
% % network_int1.eigenvector_cent(network_int1.eigenvector_cent > 0) = 1;
% cfg = [];
% cfg.filetype  = 'nifti';
% cfg.parameter = gtm;
% % cfg.filename  = './output';
% cfg.filename  = fullfile(outputdir_net,'surf.nii');
% ft_volumewrite(cfg, network_int1)
%
% restoredefaultpath
% addpath(genpath([cd_org,'/functions']));
% % close all
% addpath(connpath);
% addpath(genpath(spm_path))
% h = get(0, 'Children');
% if isempty(findobj(h,'tag','CONN functional connectivity toolbox'))
%     conn
% end
%
% filenameVOL = [];
% FSfolder = fullfile(connpath,'/utils/surf');
% sphplots = [];
% connplots = [];
% facealpha = 1;
% position = [-1 0  0];
% conn_mesh_display(fullfile(outputdir_net,'surf.nii'),[],FSfolder);

%% revert to the newer ft!
% restoredefaultpath
% addpath(genpath(cfg_main.p.ft_path));
% addpath(genpath(cfg_main.p.hcp_path));
% addpath(genpath([cfg_main.p.cd_org,'/functions']));
% addpath(genpath([cfg_main.p.cd_org,'/Data_file']));

%% parcellation - aal (132 rois)

% network_diff_lcmv.eigenvector_cent = (network_diff_lcmv.eigenvector_cent)./max(network_diff_lcmv.eigenvector_cent);
% atlas = ft_read_atlas('/data/MEG/Vahab/Github/MCW-MEGlab/tools/ft_packages/fieldtrip_20190419/template/atlas/brainnetome/BNA_MPM_thr25_1.25mm.nii');
% [~, data_intpar, coor] = vy_parcellate(network_diff_lcmv, cfg_main.atlas, gtm); % this part works with ft 2018 only!
% data_intpar.eigenvector_centdimord = 'chan';
%

%% Conn
% conn_par.conn_thre = 0.95;
% conn_ratio = vy_connvis(source_conn_pst,source_conn_bsl,conn_par, individual_headmodel, network_diff_lcmv);
% view([156,47]);

%
%% ROI summary
% [ROI, ROI_sel] = vy_ROI_report(data_intpar,.7, coor, gtm);
% disp(ROI_sel)
% savepath = fullfile(outputdir_net,['n_ROIs_',subj]);
% hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);


end