function new = do_seed_inspection(cfg_main)


mridir = cfg_main.Datalog.mridir;
subj = cfg_main.Datalog.subj;
mri_realigned = cfg_main.anat.mri_realigned;
sens = cfg_main.anat.sens;
individual_grid = cfg_main.anat.individual_grid;
individual_headmodel = cfg_main.anat.individual_headmodel;
headshape = cfg_main.anat.headshape;

%%

T1nii = fullfile(mridir,[subj(7:end), 's01_T1w.nii']); T1 = ft_read_mri(T1nii);
tran = cfg_main.anat.brain.transform;

funparam = T1.anatomy;

spherenii = fullfile(mridir,[subj(7:end), 's01_5mm_sphere.nii']);
sphere = ft_read_mri(spherenii);

A = sphere.anatomy;
[~, max_idx] = max(A(:));
[max_x, max_y, max_z] = ind2sub(size(A), max_idx);
% ft_sourceplot([], T1);
% hold on, plot3(max_x, max_y, max_z, 'm.','MarkerSize',20);


new = apply_trf(tran, [max_x, max_y, max_z]);


if cfg_main.plot == 1
    
    figure
    ft_plot_ortho(funparam, 'transform', tran, 'style', 'intersect');
    axis vis3d
    view ([-100 0])
    hold on, plot3(new(1), new(2), new(3), 'm.','MarkerSize',80);
    
    figure;
    ft_plot_mesh(individual_headmodel.bnd, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
    hold on;
    ft_plot_headshape(headshape);
    ft_plot_mesh(individual_grid.pos(individual_grid.inside, :));
    view ([-100 0])
    hold on, plot3(new(1), new(2), new(3), 'm.','MarkerSize',80);
    
    %%
    ft_determine_coordsys(mri_realigned, 'interactive', false); title('individual_mri')
    hold on
    view ([-100 0])
    ft_plot_sens(sens)
    ft_plot_mesh(individual_grid.pos(individual_grid.inside, :));
    hold on, plot3(new(1), new(2), new(3), 'm.','MarkerSize',80);
end

