function anat = do_anatomy(cfg_main)

AnatDir = cfg_main.outputmridir;
cd(AnatDir)

if exist(fullfile(cfg_main.outputmridir,['anat_',cfg_main.subj,'.mat']), 'file') == 2
    load(fullfile(cfg_main.outputmridir,['anat_',cfg_main.subj,'.mat']));
    load(fullfile(cfg_main.outputmridir,['mesh8mm_',cfg_main.subj,'.mat']));
    load(fullfile(cfg_main.outputmridir,['mesh10mm_',cfg_main.subj,'.mat']));
    
    anat = [];
    anat.mri_realigned = mri_realigned;
    anat.individual_headmodel = individual_headmodel;
    anat.headshape = headshape;
    anat.brain = brain;
    anat.individual_grid_8mm = individual_grid_8mm;
    anat.individual_grid_10mm = individual_grid_10mm;
else
    if exist(cfg_main.mripfile,'file')== 2
        
        %%
        individual_org = ft_read_mri(cfg_main.mripfile);
        individual_org = ft_convert_units(individual_org, 'mm');
        
        fid = cfg_main.fid;
        
        %% spm8 bias correction
        nii_filename_final = 'T1.nii';
        cfg                 = [];
        cfg.filename        = nii_filename_final;
        cfg.filetype        = 'nifti';
        cfg.parameter       = 'anatomy';
        ft_sourceplot([], individual_org);
        
        disp('Yes = 1, No = 2')
        bc = input('bias correction?');
        if bc ==1 %Bias correction
            ft_write_mri('BC.nii', individual_org, 'dataformat', 'nifti');
            addpath(fullfile(cfg_main.all_path.ft_path,'/external/spm8'));
            biasfield = spm_bias_estimate('BC.nii');
            spm_bias_apply('BC.nii', biasfield);
            MRI_BC = ft_read_mri('mBC.nii');
            ft_sourceplot([], MRI_BC);
            individual_org = MRI_BC;
        end
        
        disp('Yes = 1, No = 2')
        reslice_ask = input('volumereslice?');
        if reslice_ask == 1
            cfg                 = [];
            cfg.resolution      = 1;
            cfg.dim             = [256 256 256];
            individual_org_resliced2        = ft_volumereslice(cfg, individual_org);
            ft_sourceplot([], individual_org_resliced2);
        else
            individual_org_resliced2 = individual_org;
        end
        
        disp('Yes = 1, No = 2')
        slice_ask = input('is this looking acceptable?');
        if slice_ask ~= 1
            cfg             = [];
            cfg.interactive = 'yes';
            cfg.coordsys    = 'spm';
            individual_mri_spm   = ft_volumerealign(cfg, individual_org);
            cfg                 = [];
            cfg.resolution      = 1;
            cfg.dim             = [256 256 256];
            individual_org_resliced2        = ft_volumereslice(cfg, individual_mri_spm);
            ft_sourceplot([], individual_org_resliced2);
        end
        close all,
        
        %% coregister to anterior commissure based RAS space
        fprintf('Please identify the Anterior Commissure, Posterior Commissure, a point on the positive Z and X axes, and a point on the right part of the head\n');
        cfg             = [];
        cfg.interactive = 'yes';
        cfg.coordsys    = 'spm';
        individual_mri_spm   = ft_volumerealign(cfg, individual_org_resliced2);
        
        %%
        headshape = ft_read_headshape(cfg_main.hsfile);
        headshape = ft_convert_units(headshape, 'mm');
        %
        %%
        cfg = [];
        cfg.method = 'headshape';
        cfg.headshape.interactive = 'no';
        cfg.headshape.icp = 'yes';
        cfg.headshape.headshape = headshape;
        cfg.coordsys = 'neuromag';
        cfg.spmversion     = 'spm12';
        mri_realigned = ft_volumerealign(cfg, individual_mri_spm);
        
        cfg = [];
        cfg.method = 'headshape';
        cfg.headshape.interactive = 'no';
        cfg.headshape.icp = 'yes';
        cfg.headshape.headshape = headshape;
        cfg.coordsys = 'neuromag';
        cfg.spmversion     = 'spm12';
        mri_realigned = ft_volumerealign(cfg, mri_realigned);
        
        %% Subject Coordinate System (SCS / CTF)
        idx = strfind(fid.label,'LPA');
        for i=1:length(idx), idx2(i) = isempty(find(idx{i,1} == 1, 1)); end
        LPA_idx = idx2 == 0;
        idx = strfind(fid.label,'RPA');
        for i=1:length(idx), idx2(i) = isempty(find(idx{i,1} == 1, 1)); end
        RPA_idx = idx2 == 0;
        idx = strfind(fid.label,'Nasion');
        for i=1:length(idx), idx2(i) = isempty(find(idx{i,1} == 1, 1)); end
        NAS_idx = idx2 == 0;
        
        %%
        fid = ft_convert_units(cfg_main.fid,'mm');
        ft_determine_coordsys(mri_realigned, 'interactive', 'no')
        ft_plot_headshape(headshape);
        view([-90, 0]),
        hold on
        plot3(fid.pos(NAS_idx,1), fid.pos(NAS_idx,2), fid.pos(NAS_idx,3), 'm.','MarkerSize',80);
        plot3(fid.pos(LPA_idx,1), fid.pos(LPA_idx,2), fid.pos(LPA_idx,3), 'm.','MarkerSize',80);
        plot3(fid.pos(RPA_idx,1), fid.pos(RPA_idx,2), fid.pos(RPA_idx,3), 'm.','MarkerSize',80);
        
        %%
        disp('FT = 1, SPM = 2')
        seg_ask = input('volume segment?');
        switch seg_ask
            case 1
                %
                cfg = [];
                cfg.output = 'brain';
                cfg.spmversion = 'spm12';
                cfg.coordsys  = 'neuromag';
                brain = ft_volumesegment(cfg, individual_mri_spm);
            case 2
                ft_write_mri(fullfile(AnatDir, 'in_Seg.nii'), individual_mri_spm, 'dataformat', 'nifti');
                cfg = [];
                cfg.spm_path = cfg_main.all_path.spmpath;
                cfg.nii_file = 'in_Seg.nii';
                brain = do_spm_segmentation(cfg);
        end
        
        %%
        brain.transform = mri_realigned.transform;
        cfg = [];
        cfg.method = 'singleshell';
        cfg.spmversion = 'spm12';
        individual_headmodel = ft_prepare_headmodel(cfg, brain);
        
        %% Source model, warpping with template
        load temp_grid % low-res
        cfg                 = [];
        cfg.grid.warpmni    = 'yes';
        cfg.spmversion     = 'SPM12';
        cfg.grid.nonlinear  = 'yes';
        cfg.grid.template   = template_grid;
        cfg.mri             = mri_realigned;
        cfg.grid.unit       = 'mm';
        individual_grid_10mm     = ft_prepare_sourcemodel(cfg);
        
        %%
        load temp_grid_8mm % high-res
        cfg.grid.template   = template_grid;
        individual_grid_8mm     = ft_prepare_sourcemodel(cfg);
        
        %%
        % Get a list of all .nii files in the directory and, loop through each .nii file and delete it
        niiFiles = dir(fullfile(AnatDir, 'c*.nii'));
        for i = 1:numel(niiFiles)
            filePath = fullfile(AnatDir, niiFiles(i).name); delete(filePath); fprintf('Deleted file: %s\n', filePath);
        end
        
    end
    anat = [];
    anat.mri_realigned = mri_realigned;
    anat.individual_headmodel = individual_headmodel;
    anat.headshape = headshape;
    anat.brain = brain;
    anat.individual_grid_8mm = individual_grid_8mm;
    anat.individual_grid_10mm = individual_grid_10mm;
    
    save(fullfile(cfg_main.outputmridir,['anat_',cfg_main.subj,'.mat']), 'brain','mri_realigned','individual_headmodel','headshape');
    save(fullfile(cfg_main.outputmridir,['mesh8mm_',cfg_main.subj,'.mat']), 'individual_grid_8mm');
    save(fullfile(cfg_main.outputmridir,['mesh10mm_',cfg_main.subj,'.mat']), 'individual_grid_10mm');
    
end

%% Quick inspection
if cfg_main.plotflag == 1
    %%
    figure;
    ft_plot_mesh(individual_headmodel.bnd, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
    hold on;
    ft_plot_headshape(headshape);
    ft_plot_mesh(individual_grid_10mm.pos(individual_grid_10mm.inside, :));
    view ([0 90])
    
    %% plotting
    sens = ft_read_sens(cfg_main.hsfile); sens = ft_convert_units(sens, 'mm');
    figure; ft_plot_mesh(individual_headmodel.bnd, 'facecolor', 'cortex', 'edgecolor', 'none');camlight;
    hold on; ft_plot_sens(sens)
    ft_plot_headshape(headshape);
    
    %%
end

