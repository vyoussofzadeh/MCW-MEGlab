
nic = 3;
%%
savefile = [tag,'_pow.mat'];
if exist(savefile, 'file') == 2
    load(savefile)
else
    clear pow names data_dis k Sub_all
    pow = []; k = 1; 
    for i=1:length(files_sel)
        load(files_sel{i});
        datafile = files_sel{i}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
        Index = strfind(datafile, '/');
        Run = datafile(Index(end-2)+1:Index(end-1)-1);
        Subj  = datafile(Index(end-4)+1:Index(end-3)-1);
        if ~isempty(find(contains(sub_list_sel,Subj)==1, 1))
            tmp = source_ica.pow;
            tmp = tmp./max(tmp(:));
            pow(k,:) = nanmean(tmp(:,1:nic),2);
            disp(datafile)
            disp(['subj:',Subj])
            disp(['run:',Run])
            disp([num2str(k),'/', num2str(length(files_sel))])
            disp('============');
            names{k} = [num2str(k),'_',Subj,'_',Run];
            Sub_all{k} = Subj;
            k=k+1;
        end
    end
    %     save([tsk,'_pow.mat'],'pow', 'names','Sub_all','source_diff_dics','-v7.3');
end

%%
clear b
for i=1:length(names)
    b{i} = num2str(i);
end
disp(table([b',names']))

%% intra-subject
Sub_all_uniq = unique(Sub_all);

clear pow_suj
for i=1:length(Sub_all_uniq)
    Index = strfind(Sub_all,Sub_all_uniq(i));
    Index2 = find(~cellfun('isempty',Index)==1);
    pow_suj(i,:) = mean(pow(Index2,:),1);
end
subj_pow = [];
subj_pow.sub = Sub_all_uniq;
subj_pow.pow = pow_suj;

%%
template_mri = ft_read_mri(fullfile(allpath.ft_path,'template/anatomy','single_subj_T1.nii')); %
msk = 'pow';

%%
% load('/group/jbinder/ECP/MEG/group/RestEO/MNE_ICA/data/source_temp.mat')
load(fullfile(ECP_scriptdir,'source_temp.mat'))
s_in.pow = nanmean(subj_pow.pow,1)';
source_sel = s_in;

%%
if flag.map ==1
    addpath('/data/MEG/Vahab/Github/MCW-MEGlab/Projects/ECP/Run')
    Run_plot_MNEICA_subject
end
%%
% Run_plot_MNEICA_subject_stats

%%
% Run_parcellation_KB_subject

%%
% Run_plot_MNEICA_subject_run
% Run_plot_MNEICA_subject_run_stats

%% Outliers
% %         [ri, idx] = vy_outlier(evc);
% [ri, idx] = vy_outlier_baseline(pow,1);
% lout = round(length(idx)/4); % 1/5 of the outliers are thrown
% outliers_ID = (names(idx(1:lout)));
% disp('bad data:')
% disp(outliers_ID');
% disp(idx(1:lout)');
% Good_ID = (names(idx(end:-1:end-lout)));
% disp('Best data:')
% disp(Good_ID');
% disp(idx(end:-1:end-lout)');
% pow_sel = pow(idx(lout+1:end),:);
% % pow_sel = pow;
%
% %%
% idx1 = idx(lout+1:end);
% clear data_dis_sel names_sel
% for i=1:size(pow_sel,1)
%     data_dis_sel{i} = files_sel{idx1(i)};
%     names_sel{i} = names{idx1(i)};
%     sub_sel{i} = Sub_all{idx1(i)};
% end

%% normalize power values,
% for i=1:size(pow,1)
%     tmp = abs(pow(i,:));
%     tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:))); %
%     tmp(isnan(tmp)) = 0;
%     npow(i,:) = tmp;
% end

%% G-average
% evc_n = [];
% D_mean = source_sel;
% mpow = squeeze(mean(atanh(subj_pow.pow),1)); % fisher-score transformation
% % mpow = vy_normalize((pow));
% % mpow = mean(pow,1);
% % mpow = mean(pow,1);
% % D.(msk) = mpow';
%
% %%
% % D = source_sel;
% D_mean.pos     = template_grid.pos;
% D_mean.dim     = template_grid.dim;
% D_mean.inside  = template_grid.inside;
% D_mean.avg.(msk) = mpow';

%%
% [h,p,ci,stats] = ttest(pow);
% D.stat = stats.tstat';
% msk = 'stat';

%%
% close all
% cfg = [];
% cfg.mask = msk;
% cfg.loc = 'max';
% cfg.template = template_mri;
% cfg.savefile = [];
% cfg.volnorm     = 2; % yes: 1
% cfg.method  = 'ortho';
% cfg.nslices = 16;
% vy_source_plot(cfg, D_mean);
% suptitle(['Source map, mean']);
%
% %%
% close all
% cfg = [];
% cfg.parameter = msk;
% % cfg.interpmethod = 'sphere_avg';
% cfg.interpmethod = 'smudge';
% cfg.coordsys     = 'mni';
% source_int_indv = ft_sourceinterpolate(cfg, D_mean, template_mri);
% %                     source_int_indv.kurtosis(isnan(source_int_indv.kurtosis))=0;
%
% tmp = abs(source_int_indv.(msk));
% tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:))); %
% source_int_indv.(msk) = tmp;
% %
%
% cfg = [];
% cfg.subj = ('Source map, mean');
% cfg.mask = msk;
% cfg.thre = 0.7;
% cfg.savepath = [];
% cfg.colorbar = 2;
% cfg.saveflag = 2;
% cfg.surfinflated   = 'surface_inflated_both.mat';
% cfg.views = [90 0;0 90;-90 0;-180,-90];
% vy_mapvisualisation2(cfg, source_int_indv);

%%
% stat_group = vy_betweensub_stat(pow,DFN2);

%%
% if exist(fullfile(outputdir,'indiv'), 'file') == 0, mkdir(fullfile(outputdir,'indiv')); end
% cd(fullfile(outputdir,'indiv'))
% close all
% for i=1:size(pow,1)
%
%     if ~exist([names{i},'_',num2str(1),'.png'],'file')
%
%         D.avg.(msk) = pow(i,:);
%         cfg = [];
%         cfg.parameter = msk;
%         % cfg.interpmethod = 'sphere_avg';
%         cfg.interpmethod = 'smudge';
%         cfg.coordsys     = 'mni';
%         source_int_indv = ft_sourceinterpolate(cfg, D, template_mri);
%         %                     source_int_indv.kurtosis(isnan(source_int_indv.kurtosis))=0;
%
%         %
%         tmp = abs(source_int_indv.(msk));
%         tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:))); %
%         source_int_indv.(msk) = tmp;
%         %
%
%         cfg = [];
%         cfg.subj = names{i};
%         cfg.mask = msk;
%         cfg.thre = 0.5;
%         cfg.savepath = [];
%         cfg.colorbar = 2;
%         cfg.saveflag = 1;
%         cfg.surfinflated   = 'surface_inflated_both.mat';
%         cfg.views = [90 0;0 90;-90 0];
%         vy_mapvisualisation2(cfg, source_int_indv);
%
%         disp([num2str(i),'/', num2str(size(pow,1))])
%         disp(names{i})
%         pause(2)
%         close all
%     end
% end
% cd ..

%%
% close all
% savefig = fullfile(outputdir,[tag,'group']);
%
% D_mean.pos     = template_grid.pos;
% D_mean.dim     = template_grid.dim;
% D_mean.inside  = template_grid.inside;
% tmp = abs(mpow);
% tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:))); %
% mpow = tmp;
% D_mean.avg.(msk) = mpow';
% % D.avg.(msk)(isnan(D.avg.(msk)))=0;
%
% cfg = [];
% cfg.mask = msk;
% % cfg.loc = 'min';
% cfg.template = template_mri;
% cfg.savefile = [];
% cfg.volnorm     = 2; % yes: 1
% cfg.method  = 'ortho';
% cfg.nslices = 16;
% D1 = vy_source_plot(cfg, D_mean);
% set(gcf,'name','group','numbertitle','off')

%%
% savepath = fullfile(outputdir,'groupave.mat');
% cd(outputdir)
% save(savepath, 'D_mean','pow', '-v7.3');

%%
% clear savepath
% savepath{1} = fullfile(outputdir,[tag,'dics_group_2']);
% savepath{2} = fullfile(outputdir,[tag,'dics_group_3']);
%
% cfg = [];
% cfg.subj = 'group';
% cfg.mask = msk;
% cfg.thre = 0.6;
% cfg.savepath = savepath;
% cfg.saveflag = 1;
% vy_mapvisualisation(cfg, D1);

% cfg = [];
% cfg.subj = ('Source map, mean');
% cfg.mask = msk;
% cfg.thre = 0.6;
% cfg.savepath = [];
% cfg.colorbar = 2;
% cfg.saveflag = 2;
% cfg.surfinflated   = 'surface_inflated_both.mat';
% cfg.views = [90 0;0 90;-90 0;-180,-90];
% vy_mapvisualisation2(cfg, D1);

%% save group average as nii
% close all
% savenii = fullfile(outputdir,[tag,'_groupave.nii']);
% vy_savenifti(D1,msk,savenii);
%
% %% Pa
% D_mean.time = 1;
% % D_mean.avg.(msk) = abs(D_mean.avg.(msk));
% tmp = abs(mpow);
% tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:))); %
% D_mean.avg.(msk) = tmp;
%
% [~, D_par, coor] = vy_parcellate(D_mean, atlas, msk);
% % mskdim = [msk,'dimord'];
% % D_par.(mskdim) = 'chan';
%
% savepath = fullfile(outputdir,[tag,'_group_par.mat']);
% save(savepath, 'D_par', '-v7.3');

%%
% clear savepath
% savepath{1} = fullfile(outputdir,[tag,'dics_par_group_2']);
% savepath{2} = fullfile(outputdir,[tag,'dics_par_group_3']);
%
% cfg = [];
% cfg.subj = [tag,'_group'];
% cfg.mask = msk;
% cfg.thre = 0.75;
% % cfg.savepath = savepath;
% cfg.savepath = [];
% % cfg.saveflag = 1;
% cfg.saveflag = 2;
% vy_mapvisualisation(cfg, D_par);

% %%
% savenii = fullfile(outputdir,[task,'group_par.nii']);
% vy_savenifti(D_par, msk, savenii);
%
% [ROI, ROI_sel] = vy_ROI_report(D_par,.8, coor, msk);
% disp(ROI_sel)
% savepath = fullfile(outputdir,'n_ROIs_wb_par');
% hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);

% textfile_rej = fullfile(outputdir,'ROI_sel_wb_par');
% writetable(ROI_sel,textfile_rej,'Delimiter',' ');

%% Individuals
% if exist(fullfile(outputdir,'indv2'), 'file') == 0, mkdir(fullfile(outputdir,'indiv')); end
% % D = source_diff_dics;
% for i=1:size(pow,1)
%
%     D.avg.(msk) = pow(i,:)';
%
%     cfg = [];
%     cfg.mask = msk;
% %     cfg.loc = 'min';
%     cfg.template = template_mri;
%     cfg.method   = 'ortho';
%     cfg.nslices  = 20;
%     cfg.savefile = [];
%     cfg.volnorm     = 2; % yes: 1
%     vy_source_plot(cfg, D);
%     set(gcf,'Name',names{i}) %select the name you want
%     %     print(['./', tsk,'/',names{i}],'-depsc');
%     print(['./indv2/',names{i}],'-dpng');
%     disp([num2str(i),'/', num2str(size(pow,1))])
%     disp(names{i})
%     pause(2)
%     close all
%
% end

% %% Parcellation - Individuals
% disp(names')
% clear par_meg_indv
% outputinddir = fullfile(outputdir,'indiv');
% if exist(outputinddir, 'file') == 0
%     mkdir(outputinddir);   %create the directory
% end
% 
% clear savenii
% D_par = D_mean;
% 
% for k=1:size(pow,1)
%     
%     
%     disp(num2str(k))
%     Index = strfind(files_sel{k}, '/');
%     d = names{k};
%     D_par.avg.(msk) = pow(k,:)';
%     
%     %-high-res
%     cfg = [];
%     cfg.parameter = msk;
%     cfg.interpmethod = 'sphere_avg';
%     cfg.coordsys     = 'mni';
%     D1 = ft_sourceinterpolate(cfg, D_par, template_mri);
%     savenii1{k} = fullfile(outputinddir,[d(Index(1)+2:end-3),'.nii']);
%     vy_savenifti(D1,msk,savenii1{k});
%     
%     %-
%     close all
%     vol_meg = ft_read_mri(savenii1{k});
%     [~, par_meg, ~] = vy_parcellate(vol_meg, atlas, 'anatomy');
%     %     [~, par_meg, ~] = vy_parcellate(D_par, atlas, msk);
%     par_indv(k,:) = par_meg.anatomy;
%     
%     %     vy_parcellate_plot(par_meg, coor, 'net');
%     %     close all
%     
% end
% save('./par_meg','par_indv','par_meg','names','coor','savenii1')
% load('./par_meg')
% % load('par_meg');
% 
% %% Surface mapping, individuals
% addpath(allpath.connpath);
% cd indiv/
% projthresh = 0.60;
% for i=1:length(savenii1)
%     
%     nii_vol = ft_read_mri(savenii1{i});
%     vol_thre = vy_vol_thresh(nii_vol, projthresh, 'anatomy'); % abs
%     
%     Opt = [];
%     Opt.savenii = 1; Opt.savefig = 1;
%     Opt.savename = ['thre_',names{i}];
%     Opt.view = '-row';
%     vy_surfce_vis2(vol_thre,[Opt.savename,'.nii'], Opt);
%     
%     close all
% end
% cd ..
% 
% %%
% %-parcellation check
% addpath(allpath.spm_path);
% % par_meg.anatomy = mean(par_indv,1)';
% par_meg.anatomy = par_indv(16,:)';
% vy_parcellate_plot(par_meg, coor, 'net');
% 
% %%
% % par_indv1 =  par_indv;
% % par_indv1(par_indv1>0)=0;
% data = [];
% data.value = par_indv;
% data.label = par_meg.label;
% LI_meg = vy_laterality(data);
% % print(['Lat_',tsk],'-dpng')
% 
% %% LI, all subjects,
% subj_meg = names_sel;
% subj_meg = names;
% 
% figure,bar(LI_meg);
% view([90 -90])
% L = length(LI_meg);
% set(gca,'Xtick', 1:L,'XtickLabel',subj_meg);
% % set(gca,'FontSize',10,'XTickLabelRotation',90)
% set(gcf, 'Position', [1500   500   500  500]);
% grid on
% set(gca,'color','none');
% axis on
% xlim([0,L+1])
% xlabel('Subj');
% ylabel('Laterality');
% set(gcf, 'Position', [1500   500   500  500]);
% title(tsk)
% 
% %%
% % par_crm = [-0.15,0.25;
% %             0.11, 0.21;
% %             0.24, 0.4;
% %             -0.11, -0.21;
% %             0.40, 0.51;
% %             -0.31, 0.31;
% %             0.12, 0.47;
% %             0.31, 0.25;
% %             0.28, 0.43;
% %             0.13, 0.27;
% %             -0.24, -0.15
% %             0.21, 0.31];
% %
% % save('laterality_CRM','par_crm');
% % load laterality_CRM
% figure,
% % subplot 121
% hbar = bar(LI_meg);
% view([90 -90])
% L = length(LI_meg);
% for i=1:length(idx)
%     S{i} = num2str(idx(i));
% end
% % set(gca,'Xtick',idx,'XtickLabel',S);
% xlim([0,L+1]);
% ylim([-1,1]);
% set(gca,'color','none');
% box off
% xlabel('Subj');
% ylabel('Laterality');
% set(gcf, 'Position', [800   500   1000  500]);
% title(tsk,'fontsize',16)
% % legend([{'fMRI-tValue'};{'MEG'}]);
% % set(gca,'FontName','Arial');
% set(gca,'FontName','HelveticaNeueLT Std Lt');
% ylabel('Laterality','FontSize', 16);
% xlabel('Subj','FontSize', 16);
% for i=1:L
%     S{i} = num2str(i);
% end
% set(gca,'Xtick',1:L,'XtickLabel',S);
% 
% %%
% % eeg = D_par.eigenvector_cent;
% idx = [12:2:20,80:2:88];
% idx2 = [idx, idx-1];
% 
% % par_eeg.anatomy = eeg2;
% % sourceint_pow = vy_parcellate_plot(par_eeg, coor, 'net');
% % bar_input = eeg(idx2)';
% 
% bar_input = mean(par_indv,1)'; bar_input = bar_input(idx2);
% errorbar_input = std(par_indv)'; errorbar_input = errorbar_input(idx2);
% % errorbar_input=[fmri_sd(idx2),eeg_sd(idx2)]';
% 
% clear label
% for i=1:length(idx2)
%     label{i}  = par_meg.label{idx2(i)};
% end
% 
% label = {'Frontal-Inf-Oper-R' ;
%     'Frontal-Inf-Tri-R'  ;
%     'Frontal-Inf-Orb-R'  ;
%     'Rolandic-Oper-R'    ;
%     'Supp-Motor-Area-R'  ;
%     'Heschl-R'           ;
%     'Temporal-Sup-R'     ;
%     'Temporal-Pole-Sup-R';
%     'Temporal-Mid-R'     ;
%     'Temporal-Pole-Mid-R';
%     'Frontal-Inf-Oper-L' ;
%     'Frontal-Inf-Tri-L'  ;
%     'Frontal-Inf-Orb-L'  ;
%     'Rolandic-Oper-L'    ;
%     'Supp-Motor-Area-L'  ;
%     'Heschl-L'           ;
%     'Temporal-Sup-L'     ;
%     'Temporal-Pole-Sup-L';
%     'Temporal-Mid-L'     ;
%     'Temporal-Pole-Mid-L'};
% 
% L = length(bar_input);
% figure, bar(abs(bar_input), 'BarWidth', 0.9);
% % figure
% set(gca,'Xtick', 1:L,'XtickLabel',label);
% % errorbar_groups(bar_input',errorbar_input', 'bar_names',label,'bar_width',0.75,'errorbar_width',0.5);
% set(gca,'FontSize',10,'XTickLabelRotation',90);
% % grid
% box off
% set(gca,'color','none');
% set(gcf, 'Position', [500   500   1000   400]);
% set(gca,'FontSize',10,'XTickLabelRotation',90);
% set(gca,'FontName','Arial');
% % legend(tsk);
% savepath = fullfile(outputdir,'n_ROIs');
% hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);
% 
% roi = bar_input;
% save([tsk,'_ROI.mat'], 'roi', '-v7.3');

%%
% figure,bar(LI_meg);
% view([90 -90])
% L = length(LI_meg);
% set(gca,'Xtick', 1:L,'XtickLabel',names);
% % set(gca,'FontSize',10,'XTickLabelRotation',90)
% set(gcf, 'Position', [1500   500   500  500]);
% grid on
% set(gca,'color','none');
% axis on
% xlim([0,L+1])
% xlabel('Subj');
% ylabel('Laterality');
% set(gcf, 'Position', [1500   500   500  500]);
% title(tsk)

%%
% Colors = [217,217,255;89,89,89]/256;
%
% DataArray = mean(par_indv,1)';
% DataArray = DataArray(idx2);
%
% figure
% bar_handle = bar(DataArray,'grouped','BarWidth', 1);
% set(bar_handle(1),'FaceColor',Colors(2,:))
% % set(bar_handle(2),'FaceColor',Colors(2,:))
% L = length(DataArray);
% set(gca,'Xtick', 1:L,'XtickLabel',label); set(gca,'FontSize',10,'XTickLabelRotation',45)
% % set(gca,'Xtick', 1:L,'XtickLabel',([1:10,1:10]));
%
% grid off
% box off
% set(gca,'color','none');
% % axis on
% xlim([0,L+1])
% set(gca,'FontName','HelveticaNeueLT Std Lt');
% xlabel('ROI','FontSize', 20);
% ylabel('Power ERD ','FontSize', 20);
% set(gcf, 'Position', [1000   500   1000  400]);
% % title('VGA')
% legend({tsk})


%% Correlation between meg-net and fMRI
% switch task
%     case 1
%         load('par_fmri_crm.mat');load('par_meg_crm.mat');
%     case 2
%         load('par_fmri_vga.mat');load('par_meg_vga.mat');
% end
% fmri = par_fmri.anatomy./max(par_fmri.anatomy(:));
% meg = par_meg.anatomy./max(par_meg.anatomy(:));
% par_meg.anatomy = meg;
% save('par_meg','par_meg_indv','par_meg','names_sel','coor');
%

%% Laterality index
% close all
% load laterality_CRM
% figure(1),
% subplot 121
% hbar = bar(par_crm);
% view([90 -90])
% L = length(par_crm);
% for i=1:L
%     S{i} = num2str(i);
% end
% % set(gca,'Xtick',idx,'XtickLabel',S);
% xlim([0,L+1]);
% ylim([-1,1]);
% set(gca,'color','none');
% box off
% xlabel('Subj');
% ylabel('Laterality');
% % set(gcf, 'Position', [800   500   1000  500]);
% title('Word-recognition','fontsize',16)
% % legend([{'fMRI-tValue'};{'MEG'}]);
% % set(gca,'FontName','Arial');
% set(gca,'FontName','HelveticaNeueLT Std Lt');
% ylabel('Laterality','FontSize', 16);
% xlabel('Subj','FontSize', 16);
%
% DataArray = par_crm;
% Colors = [0.4922    0.0039    0.9063;0.9922    0.8672    0.0039];
% figure(2)
% subplot 121
% [xPositions, yPositions, ~, ~] = UnivarScatter(par_crm,'Label',{'fMRI','meg'},'MarkerFaceColor',Colors);
% ylabel('Laterality','FontSize', 16);
% xlabel('Modality','FontSize', 16);
% set(gca,'color','none');
% title('Word-recognition','fontsize',16)
% disp(['fmri: ',num2str(mean(par_crm(:,1))), '+-', num2str(std(par_crm(:,1)))]);
% disp(['meg: ',num2str(mean(par_crm(:,2))), '+-', num2str(std(par_crm(:,2)))]);
% set(gca,'FontName','HelveticaNeueLT Std Lt');
% hold on
% f = [xPositions, yPositions];
% for j=1:length(f)
%     line([f(j,1),f(j,2)],[f(j,3),f(j,4)]);
% end
%
%
% for k = 1:numel(hbar)
%     set(hbar(k),'FaceColor',Colors(k,:))
% end
%
% %
% load laterality_VGA
% figure(1)
% subplot 122
% hbar = bar(par_vga);
% view([90 -90])
% L = length(par_vga);
% for i=1:length(L)
%     S{i} = num2str(i);
% end
% % set(gca,'Xtick',idx,'XtickLabel',S);
% xlim([0,L+1]);
% ylim([-1,1]);
% set(gca,'color','none');
% box off
% xlabel('Subj');
% ylabel('Laterality');
% % set(gcf, 'Position', [800   500   1000  500]);
% title('Verb-Generation','fontsize',16)
% % legend([{'fmri-tValue'};{'meg-hubs'}]);
% % set(gca,'FontName','Arial');
% set(gca,'FontName','HelveticaNeueLT Std Lt');
% disp(['fmri: ',num2str(mean(par_vga(:,1))), '+-', num2str(std(par_vga(:,1)))]);
% disp(['meg: ',num2str(mean(par_vga(:,2))), '+-', num2str(std(par_vga(:,2)))]);
% % colormap(Colors);
% ylabel('Laterality','FontSize', 16);
% xlabel('Subj','FontSize', 16);
%
%
% for k = 1:numel(hbar)
%     set(hbar(k),'FaceColor',Colors(k,:))
% end
%
% figure(2)
% subplot 122
% [xPositions, yPositions, Label, RangeCut] = UnivarScatter(par_vga,'Label',{'fmri','meg'},'MarkerFaceColor',Colors);
% ylabel('Laterality','FontSize', 16);
% set(gca,'color','none');
% set(gca,'FontName','HelveticaNeueLT Std Lt');
% title('Verb-Generation','fontsize',16)
% xlabel('Modality','FontSize', 16);
% f = [xPositions, yPositions];
% for j=1:length(f)
%     line([f(j,1),f(j,2)],[f(j,3),f(j,4)]);
% end
