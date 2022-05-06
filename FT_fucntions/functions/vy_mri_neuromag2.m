function [mri_realigned,individual_headmodel,headshape, individual_grid_8mm, individual_grid_10mm] = vy_mri_neuromag2(cfg_main)

if exist(fullfile(cfg_main.outputmridir,['anat_',cfg_main.subj,'.mat']), 'file') == 2
    load(fullfile(cfg_main.outputmridir,['anat_',cfg_main.subj,'.mat']));
    load(fullfile(cfg_main.outputmridir,['mesh8mm_',cfg_main.subj,'.mat']));
    load(fullfile(cfg_main.outputmridir,['mesh10mm_',cfg_main.subj,'.mat']));
else
    if exist(cfg_main.mripfile,'file')== 2
        
        %%
        individual_mri = ft_read_mri(cfg_main.mripfile);
        %     ft_sourceplot([], individual_mri);
        individual_mri = ft_convert_units(individual_mri, 'mm');
        
        %%
        fid = cfg_main.fid;
        
        %%
        %         fprintf('Please identify the LPA, RPA, nasion, and a point on the positive Z-axis\n');
        %         cfg = [];
        %         cfg.method = 'interactive';
        %         cfg.coordsys = 'neuromag';
        %         individual_mri = ft_volumerealign(cfg, individual_mri);
        
        %% Normalized coordinate system (NCS)
        cfg          = [];
        % cfg.method   = 'interactive';
        cfg.method = 'fiducial'; % the following voxel coords were determined interactive
        cfg.coordsys = 'neuromag';
        cfg.fiducial.ac   = fid.NCS.AC;
        cfg.fiducial.pc   = fid.NCS.PC;
        cfg.fiducial.xzpoint  = fid.NCS.IH;
        cfg.fiducial.right   = fid.SCS.RPA;
        mri_realigned     = ft_volumerealign(cfg, individual_mri);
        
        %% Subject Coordinate System (SCS / CTF)
        cfg = [];
        cfg.method = 'fiducial';
        cfg.coordsys = 'neuromag';
        cfg.fiducial.nas   = fid.SCS.NAS;
        cfg.fiducial.lpa   = fid.SCS.LPA;
        cfg.fiducial.rpa   = fid.SCS.RPA;
        cfg.spmversion     = 'spm12';
        mri_realigned = ft_volumerealign(cfg, mri_realigned);
        %%
        
        %%
        %         ft_sourceplot([], mri_realigned); hold on
        %         plot3(fid.SCS.NAS(1,1), fid.SCS.NAS(1,2), fid.SCS.NAS(1,3), 'm*');
        %         plot3(fid.SCS.LPA(1,1), fid.SCS.LPA(1,2), fid.SCS.LPA(1,3), 'm*');
        %         plot3(fid.SCS.RPA(1,1), fid.SCS.RPA(1,2), fid.SCS.RPA(1,3), 'm*');
        
        %%
        headshape = ft_read_headshape(cfg_main.hsfile);
        headshape = ft_convert_units(headshape, 'mm');
        
        %%
        cfg = [];
        cfg.method = 'headshape';
        cfg.headshape.interactive = 'no';
        cfg.headshape.icp = 'yes';
        cfg.headshape.headshape = headshape;
        cfg.coordsys = 'neuromag';
        cfg.spmversion     = 'spm12';
        mri_realigned = ft_volumerealign(cfg, mri_realigned);
        
        
        %%
%         cfg = [];
%         cfg.method = 'headshape';
%         cfg.headshape.interactive = 'no';
%         cfg.headshape.icp = 'yes';
%         cfg.headshape.headshape = headshape;
%         cfg.coordsys = 'neuromag';
%         cfg.spmversion     = 'spm12';
%         mri_realigned = ft_volumerealign(cfg, mri_realigned);
        
        %% == to check everything went right with co-reg!
        ft_determine_coordsys(mri_realigned, 'interactive', 'no')
        ft_plot_headshape(headshape);
        view([-90, 0]),
        
        %%
        cfg = [];
        cfg.output = 'brain';
        cfg.spmversion = 'spm12';
        individual_seg = ft_volumesegment(cfg, mri_realigned);
        
        %%
%         cfg = [];
%         cfg.method = 'projectmesh';
%         cfg.numvertices = 10000;
%         bnd = ft_prepare_mesh(cfg, individual_seg);
        
        %%
        individual_seg.transform = mri_realigned.transform;
        cfg = [];
        cfg.method = 'singleshell';
        cfg.spmversion = 'spm12';
        individual_headmodel = ft_prepare_headmodel(cfg, individual_seg);
    end
    
    
    %% Source model, warpping with template
    load temp_grid % low-res
    % load('standard_sourcemodel3d10mm');sourcemodel = ft_convert_units(sourcemodel, 'mm');
    cfg                 = [];
    cfg.warpmni    = 'yes';
    cfg.spmversion     = 'SPM12';
    cfg.grid.nonlinear  = 'yes';
    cfg.grid.template   = template_grid;
    % cfg.grid.template   = sourcemodel;
    cfg.mri             = mri_realigned;
    cfg.grid.unit       = 'mm';
    individual_grid_10mm     = ft_prepare_sourcemodel(cfg);
    
    %%
     
    %%
    % load('standard_sourcemodel3d8mm');sourcemodel = ft_convert_units(sourcemodel, 'mm');
    clear template_grid
    load temp_grid_8mm % high-res
    cfg.grid.template   = template_grid;
    % cfg.grid.template   = sourcemodel;
    cfg.grid.template   = template_grid;
    individual_grid_8mm     = ft_prepare_sourcemodel(cfg);
    
    %%
    save(fullfile(cfg_main.outputmridir,['anat_',cfg_main.subj,'.mat']),'individual_seg','mri_realigned','individual_headmodel','headshape');
    save(fullfile(cfg_main.outputmridir,['mesh8mm_',cfg_main.subj,'.mat']), 'individual_grid_8mm');
    save(fullfile(cfg_main.outputmridir,['mesh10mm_',cfg_main.subj,'.mat']), 'individual_grid_10mm');
    
    
end
%%
% cfg = [];
% cfg.method = 'headshape';
% cfg.headshape.interactive = 'no';
% cfg.headshape.icp = 'yes';
% cfg.headshape.headshape = headshape;
% cfg.coordsys = 'ctf';
% cfg.spmversion     = 'spm12';
% mri_realigned_ctf = ft_volumerealign(cfg, mri_realigned);

%     end

%% Quick inspection
if cfg_main.plotflag == 1
    
    %%
    figure;
    ft_plot_vol(individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
    hold on;
    ft_plot_headshape(headshape);
    ft_plot_mesh(individual_grid_10mm.pos(individual_grid_10mm.inside, :));
    view ([0 90])
    
    figure;
    ft_plot_vol(individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
    hold on;
    ft_plot_headshape(headshape);
    ft_plot_mesh(individual_grid_8mm.pos(individual_grid_8mm.inside, :));
    view ([0 90])
    
    %     figure
    %     ft_plot_vol(individual_headmodel, 'unit', 'mm');  %this is the brain shaped head model volume
    % %     ft_plot_sens(t_data.all.grad, 'unit', 'mm', 'coilsize', 10);  %this is the sensor locations
    %     ft_plot_mesh(individual_grid_10mm.pos(individual_grid_10mm.inside, :));
    %     ft_plot_ortho(mri_realigned.anatomy, 'transform', mri_realigned.transform, 'style', 'intersect');
    
    %%
    sens = ft_read_sens(cfg_main.hsfile); sens = ft_convert_units(sens, 'mm');
    figure; ft_plot_vol(individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none'); camlight;
    hold on; ft_plot_sens(sens)
    ft_plot_mesh(individual_grid_8mm.pos(individual_grid_8mm.inside, :));
    ft_plot_headshape(headshape);
    
    %%
end



%%
% create the individual grid from mni grid
%
%         temp = load('standard_sourcemodel3d8mm');
%         template_grid = ft_convert_units(template_grid, 'mm');
%         cfg = [];
%         cfg.grid.warpmni = 'yes';
%         %     cfg.grid.template = temp.sourcemodel;
%         cfg.grid.template = template_grid;
%         cfg.grid.nonlinear = 'yes';
%         cfg.mri = mri_realigned;
%         cfg.grid.unit     = 'mm';
%         individual_grid = ft_prepare_sourcemodel(cfg);

%% GM masking
%
% atlas = cfg_main.atlas;
% cfg = [];
% cfg.atlas      = atlas;
% cfg.roi        = atlas.tissuelabel;  % here you can also specify a single label, i.e. single ROI
% cfg.inputcoord = 'mni';
% mask           = ft_volumelookup(cfg, individual_grid_8mm);
%
% %% After
% individual_grid2                 = individual_grid_8mm;
% individual_grid2.inside          = false(individual_grid2.dim);
% individual_grid2.inside(mask==1) = true;
%
% % figure; hold on;
% % ft_plot_vol(individual_headmodel, 'edgecolor', 'none', 'facealpha', 0.4);
% % ft_plot_mesh(individual_grid_10mm.pos(individual_grid_10mm.inside,:));
%
%
% figure;
% ft_plot_vol(individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
% hold on;
% % ft_plot_headshape(headshape);
% ft_plot_mesh(individual_grid2.pos(individual_grid2.inside, :));
% view ([0 90])


%%
%         cfg                 = [];
%         cfg.grad            = cfg_main.megdata.grad;
%         cfg.headmodel       = individual_headmodel;
%         cfg.reducerank      = 2;
%         cfg.channel         = {'MEG'};
%         cfg.resolution = 10;   % use a 3-D grid with a 1 cm resolution
%         cfg.sourcemodel.unit       = 'mm';
%         [grid] = ft_prepare_leadfield(cfg);

%     figure;
%     ft_plot_vol(individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
%     hold on;
%     ft_plot_headshape(headshape);
%     ft_plot_mesh(grid.pos(grid.inside, :));
%     view ([0 90])