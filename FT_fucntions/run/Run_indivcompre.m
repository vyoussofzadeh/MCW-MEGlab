% clear, 
close all,
clc
cd(outputdir)

%%
PN = load('./PN/par_meg.mat');
PN2 = load('./PN/PN_pow');
PN3 = load('./PN/PN_ROI');

DFN = load('./DFN/par_meg.mat');
DFN2 = load('./DFN/DFN_pow');
DFN3 = load('./DFN/DFN_ROI');

coor = PN.coor;

%%
outputdir_idivcompr = './indivcomp';
if exist(outputdir_idivcompr, 'file') == 0, mkdir(outputdir_idivcompr); end
cd(outputdir_idivcompr)

%%
[C,ia,ib] = intersect(DFN2.Sub_all,PN2.Sub_all, 'stable');
disp(C')
a = [ia,ib];

par_DFN_sel  = DFN.par_indv(ia,:);
par_PN_sel = PN.par_indv(ib,:);

pow_DFN_sel = DFN2.pow(ia,:);
pow_PN_sel  = PN2.pow(ib,:);

L = size(pow_DFN_sel,1);

%% 
template_mri = ft_read_mri(fullfile(allpath.ft_path,'template/anatomy','single_subj_T1.nii')); %
msk = 'pow';
projthresh = 0.75;

%%
addpath(allpath.connpath);
addpath(allpath.spm_path);

%% DFN (modified)
close all

tsk = 'dfn'; pow_sel = pow_DFN_sel;
D = PN2.source_diff_dics;
for i=1:size(pow_sel,1)
    
    D.(msk) = pow_sel(i,:)+mean(pow_sel,1);
    cfg = [];
    cfg.mask = 'pow';
    cfg.loc = 'min';
    cfg.template = template_mri;
    cfg.savefile = [];
    cfg.volnorm     = 2; % yes: 1
    source_dics = vy_source_plot(cfg, D);
    
    source_dics1 = vy_vol_thresh(source_dics,projthresh,'pow'); % abs
    source_dics1.pow = -(source_dics1.pow./max(source_dics1.pow(:)));
    
    set(gcf,'Name',C{i}) %select the name you want
    savenname = [tsk,'_',num2str(i),'_',C{i}];
    vy_savenifti(source_dics1,'pow',[savenname,'.nii']);
%     print(['./', tsk,'/',names{i}],'-depsc');
%     print(savenname,'-dpng');

    %-Surf-vis
    opt.run = [];
    opt.tsk = tsk;
    opt.subj = C{i};
    opt.savedir = [];
    opt.savenii = 2;
%     opt.plot = '-mosaic';
    opt.plot = '-row';
    vy_surfce_vis([],[savenname,'.nii'], opt);

    disp([num2str(i),'/', num2str(size(pow_sel,1))])
    pause(1)
    close all
end

%% For manual saving ...
% view([-110,20])
% view([110,20])

%% PN
close all

tsk = 'pn'; pow_sel = pow_PN_sel;
D = PN2.source_diff_dics;
for i=1:size(pow_sel,1)
    D.(msk) = pow_sel(i,:)+mean(pow_sel,1);
    cfg = [];
    cfg.mask = 'pow';
    cfg.loc = 'min';
    cfg.template = template_mri;
    cfg.savefile = [];
    cfg.volnorm     = 2; % yes: 1
    source_dics = vy_source_plot(cfg, D);
    source_dics1 = vy_vol_thresh(source_dics,projthresh,'pow'); % abs
    source_dics1.pow = -(source_dics1.pow./max(source_dics1.pow(:)));
    
    set(gcf,'Name',C{i}) %select the name you want
    savenname = [tsk,'_',num2str(i),'_',C{i}];
    vy_savenifti(source_dics1,'pow',[savenname,'.nii']);
    %     print(['./', tsk,'/',names{i}],'-depsc');
    %     print(savenname,'-dpng');
    
    %-Surf-vis
    opt.run = [];
    opt.tsk = tsk;
    opt.subj = C{i};
    opt.savedir = [];
    opt.savenii = 2;
    opt.plot = '-row';
    vy_surfce_vis([],[savenname,'.nii'], opt);
    
    disp([num2str(i),'/', num2str(size(pow_sel,1))])
    pause(1)
    close all
end

%% Overlay PN-DFN
projthresh_meg = 0.9;
for i=1:L
    %     for j=1:length(subj_fmri)
    %         if strncmp(subj_fmri{j},subj_meg{i},5)==1
    
    disp(C{i})
    savenname_dfn = ['dfn_',num2str(i),'_',C{i},'.nii'];
    savenname_pn = ['pn_',num2str(i),'_',C{i},'.nii'];
    
    %---
    s_vol_dfn = ft_read_mri(savenname_dfn);
    s_vol_pn = ft_read_mri(savenname_pn);
    
%     s_vol_dfn = vy_vol_thresh(s_vol_dfn,projthresh_meg,'anatomy'); % abs
%     s_vol_pn = vy_vol_thresh(s_vol_pn,projthresh_meg,'anatomy'); % abs
    
%     s_vol_pn.anatomy(isnan(s_vol_fmri.anatomy(:)))=0;
%     
%     
%     s_vol_meg = s_vol_fmri;
%     s_vol_meg.anatomy = squeeze(vol_meg(i,:,:,:));
%     
%     a = (s_vol_fmri1.anatomy); a = a./max(a(:));
%     b = (s_vol_meg.anatomy); b = b./max(b(:));
%     s_vol_meg = vy_vol_thresh(s_vol_meg,projthresh_meg); % abs
    
    savenname_comp = ['comb_',num2str(i),'_',C{i},'.nii'];
    spm_imcalc({savenname_dfn, savenname_pn},savenname_comp, 'i1.*(i1<0)-i2.*(i2<0)');
    
    %-Surf-vis
    opt.run = [];
    opt.tsk = 'comb';
    opt.subj = ['S',num2str(i),'_',C{i}];
    opt.savedir = [];
    opt.savenii = 2;
    opt.plot = '-row';
    vy_surfce_vis([],savenname_comp, opt);
    %     masked{k} = subj_fmri{j}; k=k+1;
    
    %         end
    %     end
    close all
end

%%
% cfg = [];
% cfg.mask = 'stat';
% cfg.loc  = 'min';
% cfg.template = template_mri;
% cfg.savefile = savefig;
% cfg.volnorm  = 2; % yes: 1
% source_stat_dfn = vy_source_plot(cfg, stat_group_DFN);
% 
% savenii_dnf = 'groupave_dfn.nii';
% vy_savenifti(source_stat_dfn,'stat',savenii_dnf);
% 
% %% Individuals
% if exist(fullfile(outputdir,'indiv'), 'file') == 0, mkdir(fullfile(outputdir,'indiv')); end
% D = source_diff_dics;
% for i=70:size(pow,1)
%     D.(msk) = pow(i,:)';
%     cfg = [];
%     cfg.mask = 'pow';
%     cfg.loc = 'min';
%     cfg.template = template_mri;
%     cfg.savefile = [];
%     cfg.volnorm     = 2; % yes: 1
%     vy_source_plot(cfg, D);
%     set(gcf,'Name',names{i}) %select the name you want
% %     print(['./', tsk,'/',names{i}],'-depsc');
%     print(['./', tsk,'/',names{i}],'-dpng');
%     disp([num2str(i),'/', num2str(size(pow,1))])
%     disp(names{i})
%     pause(2)
%     close all
% end

%%
