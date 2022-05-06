
cd(outputdir)

cd(outputdir)
outputdir_atlas = './atlas_ROIs';
if exist(outputdir_atlas, 'file') == 0, mkdir(outputdir_atlas); end
cd(outputdir_atlas)

%%
addpath(allpath.connpath);
addpath(allpath.spm_path);

%%
% template_mri = ft_read_mri(fullfile(allpath.ft_path,'template/anatomy','single_subj_T1.nii')); %
AAL_atlas = fullfile(allpath.connpath,'utils/otherrois/aal_ft.nii');

%%
% PN2 = load(fullfile(outputdir,'PN/PN_pow'));
% % savenii = [tsk,'_group_par_nii.nii'];
% 
% %-
% tsk = 'PN'; msk = 'pow';
% [~, D_par_PN_nii, coor] = vy_parcellate(PN2.source_diff_dics, atlas, msk);

%%
idx = [];
idx.central = [2,20]; 
idx.frontal = 4:2:18; 
idx.subcor = [22:2:48,78]; 
idx.Occ = 50:2:56; 
idx.pari= 58:2:70; 
idx.temp = 80:2:90;

% idx = [idx_central,idx_frontal,idx_subcor,idx_Occ,idx_pari,idx_temp];
idx1         = [idx.frontal,idx.temp,idx.pari];
% idx1         = [idx.frontal];
% idx1         = [idx.frontal, idx.temp];
% idx1 = 80;
idx2        = [idx1-1, idx1];
Run_aal_labels

clear label
for i=1:length(idx2)
    label{i}  = aal_label1{idx2(i)};
end
disp(label')


%%
idx = [];
idx.central = [2,20]; 
idx.frontal = 4:2:18;
idx.frontal1 = [12,14,16,18]; 
idx.subcor = [22:2:48,78]; 
idx.Occ = 50:2:56; 
idx.pari= 58:2:70;
idx.pari1= 58:2:66; 
idx.temp = 80:2:90;
idx.temp1 = 82:2:88; % no Heschl's G


% idx = [idx_central,idx_frontal,idx_subcor,idx_Occ,idx_pari,idx_temp];
idx1         = [idx.frontal1,idx.temp1,idx.pari1];
% idx1         = [idx.frontal1,idx.temp1,];
% idx1         = [idx.temp1];
% idx1         = [idx.frontal1];
% idx1         = [idx.frontal, idx.temp];
% idx1 = 80;
idx2        = [idx1-1, idx1];
Run_aal_labels

clear label
for i=1:length(idx2)
    label{i}  = aal_label1{idx2(i)};
end
disp(label')

%%
% xmin=-1;
% xmax=1;
% n=1;
% x=xmin+rand(1,n)*(xmax-xmin);
% 
% m = zeros(116,1);
% for i=1:length(idx2)
%     m(idx2(i)) = rand(1);
% %     m(idx2(i)) = 1;
%     m(idx2(i)) = xmin+rand(1,n)*(xmax-xmin);
% end
% D_par_PN_nii.pow = m;
% D_par_PN_nii.([msk,'dimord']) = 'chan';
% 
% savenii = 'ROIs_mask.nii';
% vy_savenifti(D_par_PN_nii, msk, savenii);
% 
% % mask = ft_read_mri(savenii);
% % projthresh = 0.70;
% % mask1 = vy_vol_thresh(mask, projthresh, 'anatomy'); % abs
% 
% % Opt = [];
% % Opt.savenii = 1; Opt.savefig = 0;
% % Opt.savename = 'ROIs_mask1';
% % Opt.view = '-mosaic';
% % vy_surfce_vis2(mask1,[Opt.savename,'.nii'], Opt);
% 
% %-
% % input_nii = load(savenii);
% %
% % %%
% % savenii = 'DFN_group_par_thre.nii';
% % vy_savenifti(D_par_PN_nii, msk, savenii);
% filenameSURF=vy_conn_mesh_display_expand(savenii);
% conn_mesh_display(savenii);

%%

filenameSURF = vy_conn_mesh_display_expand(AAL_atlas, idx2);
% close all
filenameSURF = vy_conn_mesh_display_expand(AAL_atlas, idx2(19));
conn_mesh_display(filenameSURF);
view([90,0])
% view([-90,0])
% view([-110,20])

%%
test = ft_read_mri('/data/MEG/Vahab/Github/MCW-MEGlab/tools/Conn/conn/rois/atlas.nii');
atlas = ft_read_atlas('/data/MEG/Vahab/Github/MCW-MEGlab/tools/Conn/conn/rois/atlas.nii');




