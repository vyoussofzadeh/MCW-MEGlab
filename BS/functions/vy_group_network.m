%%
clear evc conn names data_dis k
k = 1;
for i=1:length(files_sel)
    load(files_sel{i});
    datafile = files_sel{i}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
    Index = strfind(datafile, '/');
    Date = datafile(Index(5)+1:Index(6)-1);
    Subj  = datafile(Index(6)+1:Index(7)-1);
    pow(i,:) = source_diff_dics.pow;
    disp(datafile)
    disp(['subj:',Subj])
    disp(['Date:',Date])
    disp('============');
    names{i} = [num2str(i),'_',Subj,'_',Date];
    %     Sub_all{k} = Subj;
end
outputdir = fullfile(DestDirectory,'group_dics');
if exist(outputdir, 'file') == 0
    mkdir(outputdir);   %create the directory
end
cd(outputdir)

%%

% clear evc conn names data_dis k
% k = 1;
% for i=1:length(files_sel)
%     kk=1;
%     disp(files_sel{i})
% %     d1 = dir(fullfile(DestDirectory,files_sel{i})); files = {d1.name}'; files_sel1 = files(3:end,1);
%     for j=1:length(files_sel1)
%         d2 = fullfile(DestDirectory,files_sel{i},files_sel1{j},'network');
%         files = dir(fullfile(d2,'network_evc_C-*.mat'));
%         if isempty(files) == 0
%             load(fullfile(files.folder,files.name));
%             evc(k,:) = network_diff_lcmv.eigenvector_cent;
%             data_dis{k} = fullfile(files.folder,files.name);
%             Index = strfind(data_dis{1}, '.mat');
%             names{k} = [files_sel{i},'_',num2str(kk)];
%             sub{k} = files_sel{i};
%             disp(kk)
%             k=k+1; kk=kk+1;
%         end
%     end
% end
% outputdir = fullfile(DestDirectory,'group_network');
% if exist(outputdir, 'file') == 0
%     mkdir(outputdir);   %create the directory
% end
% cd(outputdir)

%%
clear b
for i=1:length(names)
    b{i} = num2str(i);
end
disp(table([b',names']))

%%
% evc_sel = evc;
% data_dis_sel = data_dis;
% names_sel = names;
% sub_sel = sub;

%% Individuals
template_mri = ft_read_mri(fullfile(allpath.ft_path,'template/anatomy','single_subj_T1.nii')); %
msk = 'pow';
% D = source_diff_dics;
% for i=1:size(pow,1)
%     D.(msk) = pow(i,:)';
%     cfg = [];
%     cfg.mask = 'pow';
%     cfg.loc = 'min';
%     cfg.template = template_mri;
%     cfg.savefile = [];
%     cfg.volnorm     = 2; % yes: 1
%     vy_source_plot(cfg, D);
%     
%     set(gcf,'Name',names{i}) %select the name you want
%     pause(1)
% end

%% G-average
evc_n = [];
D = source_diff_dics;
evcn = squeeze(mean(pow,1)); % zero-mean, divide by std, average
% evcn = vy_normalize(pow);

D.(msk) = evcn';

cfg = [];
cfg.mask = 'pow';
cfg.loc = 'min';
cfg.template = template_mri;
cfg.savefile = [];
cfg.volnorm     = 2; % yes: 1
dics_pow = vy_source_plot(cfg, D);

%%
savepath = fullfile(outputdir,'groupave.mat');
cd(outputdir)
save(savepath, 'D','pow', '-v7.3');

%% Source vis - quick inspection
% cfg              = [];
% cfg.method       = 'ortho';
% cfg.funparameter = msk;
% % figure
% ft_sourceplot(cfg,D);
% hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);


cfg = [];
cfg.subj = 'group';
cfg.mask = 'pow';
cfg.thre = 0.6;
cfg.savepath = savepath;
vy_mapvisualisation(cfg, D);
% vy_mapvisualisation(D,msk,0.7, []);

%%
% D.eigenvector_cent(D.eigenvector_cent < 0.7.*max(D.eigenvector_cent)) = NaN; % negative effects
% D.eigenvector_centdimord = 'chan';

%% Source vis - template, quick and dirty visualisation
% D.eigenvector_cent = 100.*D.eigenvector_cent;
% gtm = 'eigenvector_cent';
% clear savepath
% 
% savepath{1} = fullfile(outputdir,'n_par1_t_2_wb');
% savepath{2} = fullfile(outputdir,'n_par2_t_3_wb');
savenii = fullfile(outputdir,'groupave.nii');
% 
% param = [];
% param.mask = 'eigenvector_cent';
% param.loc = 'max';
% D1 = vy_source_plot(D,template_mri,param,2);
% vy_mapvisualisation(D1,'eigenvector_cent',0.6, savepath);
vy_savenifti(D,msk,savenii);

%% more detailled visualisation
[~, D_par, coor] = vy_parcellate(D, atlas, msk);
D_par.powdimord = 'chan';

savepath = fullfile(outputdir,'group_par.mat');
save(savepath, 'D_par', '-v7.3');

% clear savepath
% savepath{1} = fullfile(outputdir,'n_par1_t_2_wb_par');
% savepath{2} = fullfile(outputdir,'n_par2_t_3_wb_par');
vy_mapvisualisation(D_par,msk,0.6, []);

savenii = fullfile(outputdir,'group_par.nii');
vy_savenifti(D_par, msk, savenii);

[ROI, ROI_sel] = vy_ROI_report(D_par,.8, coor, msk);
% disp(ROI_sel)
% savepath = fullfile(outputdir,'n_ROIs_wb_par');
% hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);

% textfile_rej = fullfile(outputdir,'ROI_sel_wb_par');
% writetable(ROI_sel,textfile_rej,'Delimiter',' ');

%% Parcellation - Individuals
% clear par_meg_indv
% outputinddir = fullfile(outputdir,'indiv_n_wb');
% if exist(outputinddir, 'file') == 0
%     mkdir(outputinddir);   %create the directory
% end
% clear savenii
% for k=1:size(evc_sel,1)
%     
%     Index = strfind(data_dis_sel{k}, '\');
%     d = data_dis_sel{k};
%     D.eigenvector_cent = evc_sel(k,:)';
%     
%     %-high-res
%     cfg = [];
%     cfg.parameter = gtm;
%     cfg.interpmethod = 'sphere_avg';
%     cfg.coordsys     = 'mni';
%     D1 = ft_sourceinterpolate(cfg, D, template_mri);
%     savenii1{k} = fullfile(outputinddir,[d(Index(end)+1:end-11),names_sel{k},'.nii']);
%     vy_savenifti(D1,gtm,savenii1{k});
%     
%     %-
%     close all
%     vol_meg = ft_read_mri(savenii1{k});
%     [~, par_meg, ~] = vy_parcellate(vol_meg, atlas,'anatomy');
%     par_meg_indv(k,:) = par_meg.anatomy;
%     
% end
% save('par_meg','par_meg_indv','par_meg','names_sel','coor')
% % load('par_meg');
% 
% %-parcellation check
% addpath(spm_path);
% par_meg.anatomy = mean(par_meg_indv,1)';
% vy_parcellate_plot(par_meg, coor, 'net');

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
