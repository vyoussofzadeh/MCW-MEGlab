
% stimcoor_1 = '/data/MEG/Research/logopenicppa/MRI/idsc9005/stimulation_rois/002s01_T1w_RAI.1D';
% stimcoor_2 = '/data/MEG/Research/logopenicppa/MRI/idsc9005/stimulation_rois/004s01_T1w_RAI.1D';
% stimcoor_3 = '/data/MEG/Research/logopenicppa/MRI/idsc9005/stimulation_rois/005s01_T1w_RAI.1D';
% stimcoor_4 = '/data/MEG/Research/logopenicppa/MRI/idsc9005/stimulation_rois/006s01_T1w_RAI.1D';







% s01T1wRAI_1 = load(stimcoor_1);
% s01T1wRAI_2 = load(stimcoor_2);
% s01T1wRAI_3 = load(stimcoor_3);
% s01T1wRAI_4 = load(stimcoor_4);

sphere = ft_read_mri('/data/MEG/Research/logopenicppa/MRI/idsc9005/stimulation_rois/002s01_5mm_sphere.nii');
A = sphere.anatomy; [~, max_idx] = max(A(:)); [max_x, max_y, max_z] = ind2sub(size(A), max_idx);

anat_1 = load('/data/MEG/Research/logopenicppa/ft_process/32037_002/anat/anat_32037_002.mat');
tran = anat_1.mri_realigned.transform;
% src1 = ft_transform_geometry(tran, src);

new_1 = apply_trf(tran, [max_x, max_y, max_z]);



src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Atlas/cortex_pial_low.fs';
src = ft_read_headshape(src_fname);
vertexcolor = zeros(size(src.pos,1), 3);

% mri_temp = ft_read_mri('single_subj_T1.nii');

individual_headmodel = anat_1.individual_headmodel;
headshape = anat_1.headshape;
individual_grid = anat_1.individual_grid;

figure;
ft_plot_mesh(individual_headmodel.bnd, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
hold on;
ft_plot_headshape(headshape);
% ft_plot_mesh(individual_grid.pos(individual_grid.inside, :));
view ([-100 0])
hold on, plot3(new_1(1), new_1(2), new_1(3), 'm.','MarkerSize',80);


% close all
% figure, ft_plot_mesh(src, 'vertexcolor', vertexcolor);
% % ft_plot_mesh(src1, 'maskstyle', 'opacity', 'facecolor', 'black', 'facealpha', 0.25, 'edgecolor', 'red',   'edgeopacity', 0.5);
% 
% alpha(0.1)
% coor = new_1; hold on, plot3(-coor(1), coor(2), coor(3), 'm.','MarkerSize',80);
% coor = s01T1wRAI_2; hold on, plot3(-coor(1), coor(2), coor(3), 'm.','MarkerSize',80);
% coor = s01T1wRAI_3; hold on, plot3(-coor(1), coor(2), coor(3), 'm.','MarkerSize',80);
% coor = s01T1wRAI_4; hold on, plot3(-coor(1), coor(2), coor(3), 'm.','MarkerSize',80);
% 
% 
% sourcemodel = ft_read_headshape('cortex_20484.surf.gii')
% vertexcolor = zeros(size(sourcemodel.pos,1), 3);
% 
% close all
% figure, ft_plot_mesh(sourcemodel, 'vertexcolor', vertexcolor);
% alpha(0.1)
% 
% 
% cfg = [];
% cfg.grad = grad;
% cfg.headmodel = headmodel;
% cfg.grid = grid;
% cfg.channel = {'MEG'};
% lf = ft_prepare_leadfield(cfg);
