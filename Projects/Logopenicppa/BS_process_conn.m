%% Logopenicppa dataset, Medical College of Wisconsin

% Script: BS Process (seed-based conn analysis)
% Project: Logopenicppa_rest
% Writtern by: Vahab YoussofZadeh
% Update: 04/07/2023

clear; clc, close('all'); warning off,

%% Paths
restoredefaultpath
script_path = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/Logopenicppa';
cd(script_path)

% indir = '/data/MEG/Research/logopenicppa/raw';
% outdir = '/data/MEG/Research/logopenicppa/BS_process';
indir = '/data/MEG/Research/logopenicppa';
outdir = fullfile(indir,'ft_process');
rawdatadir = fullfile(indir,'raw');

ft_path ='/opt/matlab_toolboxes/ft_packages/Stable_version/fieldtrip-master';
addpath(ft_path); ft_defaults

%%
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/brewermap')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/functions_new')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/Colormaps-from-MatPlotLib2.0')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/Miscellaneous')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/Logopenicppa/functions')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/Data_file')

%%
all_path = [];
all_path.ft18 = fullfile('/opt/matlab_toolboxes/ft_packages/fieldtrip_041718'); % needed for IC plotting
all_path.ft_path = ft_path;
all_path.script_path = script_path;

% %- Adding path
% cfg_init = [];
% cfg_init.path_tools = '/data/MEG/Vahab/Github/MCW_MEGlab/tools';
% [allpath, atlas] = do_init(cfg_init);

%- atlas
atlas_path = fullfile(ft_path,'template','atlas');
atlas = ft_read_atlas(fullfile(atlas_path,'aal/ROI_MNI_V4.nii'));

%% Loading up raw data
clc
tag = 'Rest';
cfg = [];
cfg.indir = rawdatadir;
cfg.tag = tag;
datafile = do_datalookup(cfg);

datafile.datafile_fif
datafile.run
datafile.task
datafile.task_run

%%
cfg = []; cfg.layout = 'neuromag306mag.lay'; lay = ft_prepare_layout(cfg);

%%
mridir = '/data/MEG/Research/logopenicppa/MRI/idsc9005/stimulation_rois';

%%
for i = 1:length(datafile.task_run)
    
    datafile_sel = datafile.datafile_fif{i}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
    Index = strfind(datafile_sel, '/');
    subj = datafile_sel(Index(6)+1:Index(7)-3);
    Index = strfind(datafile_sel, 'run');
    run  = datafile_sel(Index+3);
    
    mripfile = fullfile(mridir,[subj(7:end-2), 's01_T1w.nii']);
    if exist(mripfile,'file')== 2
        
        disp(datafile_sel)
        disp(['subj:',subj])
        disp(['Run:',run])
        
        %-elec/grad
        %         sens = ft_read_sens(datafile_sel);
        %         sens = ft_convert_units(sens,'mm');
        %         disp('============');
        
        subdir = fullfile(outdir,subj, tag, run);
        if exist(subdir, 'file') == 0
            mkdir(subdir);   %create the directory
        end
        cd(subdir)
        disp(['outputdir:',subdir])
        disp('============');
        
        %% Updating datalog
        Datalog = [];
        Datalog.subj = subj;
        Datalog.run = run;
        Datalog.datafile = datafile;
        Datalog.subdir = subdir;
        Datalog.outdir = outdir;
        
        %% Preprocesssing
        flag = [];
        flag.preproces.filtering = 1;
        flag.preproces.artifact = 1;
        flag.preproces.ica = 1;
        flag.ic_selection = 1;
        
        cfg = [];
        cfg.datalog = Datalog;
        cfg.flag = flag;
        cfg.lay = lay;
        cfg.all_path = all_path;
        cln_data = do_preprocess_rest(cfg, datafile_sel);
        addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/functions_new')
        
        %% Preprocessing, notch filtering
        cfg = [];
        cfg.subj = subj;
        cfg.pad = 4;
        fn_data = do_notch(cfg, cln_data);
        
        %%
        cfg = [];
        cfg.channel = 'MEG';
        cfg.covariance = 'yes';
        cov_matrix = ft_timelockanalysis(cfg, fn_data);
        
        %%
        outputmridir = fullfile(outdir,subj,'anat'); % output dir
        if exist(outputmridir, 'file') == 0, mkdir(outputmridir); end
        cd(outputmridir)
        
        %         mridir = fullfile(indir,subj,'Anatomy/mri');
        %         fidfile = fullfile(indir,subj,'Anatomy/bem/Anatomy-fiducials.fif');
        %
        %         clear fid
        %         mripfile = fullfile(mridir,'T1.mgz');
        %         if exist(outputmridir, 'file') == 0, mkdir(outputmridir); end
        
        headshape = ft_read_headshape(datafile_sel);
        
        cfg = [];
        cfg.megdata = cov_matrix.grad;
        cfg.mripfile = mripfile;
        cfg.hsfile = datafile_sel; % headshape;
        cfg.fid = headshape.fid;
        cfg.outputmridir = outputmridir;
        cfg.subj = subj;
        cfg.plotflag = 1;
        cfg.atlas = atlas;
        [mri_realigned,individual_headmodel,headshape, individual_grid_8mm, individual_grid_10mm] = do_anatomy(cfg);
        cd(outd.sub)
        
        figure,
        ft_plot_mesh(individual_headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');
        figure; ft_plot_mesh(individual_headmodel.bnd, 'facecolor', 'cortex', 'edgecolor', 'none');
        alpha 0.5; camlight;
        
        %% Choosing mesh
        % Volumetric-based analysis
        flag.anatomy_check = 2;
        flag.meshgrid_sel = 1;
        choose_grid = 2;
        switch choose_grid
            % if flag.warping == 1
            case 1
                switch flag.meshgrid_sel
                    case 1
                        individual_grid = individual_grid_10mm;
                    case 2
                        individual_grid = individual_grid_8mm;
                end
            case 2
                switch flag.meshgrid_sel
                    case 1
                        meshtag = 'lowres';
                        load temp_grid % low-res
                        template_grid = ft_convert_units(template_grid, 'mm');
                        individual_grid = individual_grid_10mm;
                    case 2
                        meshtag = 'highres';
                        load temp_grid_8mm % high-res
                        individual_grid = individual_grid_8mm;
                end
        end
        
        %% Anatomoy check!
        saveflag = 2;
        if flag.anatomy_check == 1
            cfg = [];
            cfg.saveflag = [];
            cfg.headmodel = individual_headmodel;
            cfg.leadfield = individual_grid;
            cfg.mri_realigned  = mri_realigned;
            cfg.headshape = headshape;
            cfg.outputmridir = outputmridir;
            cfg.mtd = 'vol';
            do_mri_inspection(cfg, cln_data);
        end
        
        %%
        anat = [];
        anat.individual_grid = individual_grid;
        anat.individual_headmodel = individual_headmodel;
        anat.template_grid = template_grid;
        anat.atlas = atlas;
        
        cfg = [];
        cfg.cov_matrix = cov_matrix;
        cfg.anat = anat;
        do_wPLIconn(cfg, cln_data)
        
        %%
        close all
        funparam = mri_realigned.anatomy;
        data  = ft_checkdata(mri_realigned, 'hasunit', 'yes');
        [corner_vox, corner_head] = cornerpoints(data.dim, data.transform);
        diagonal_head = norm(range(corner_head));
        diagonal_vox  = norm(range(corner_vox));
        resolution    = diagonal_head/diagonal_vox; % this is in units of "data.unit"
        clim = [0,1];
        
        stimcoor = fullfile(mridir,[subj(7:end), 's01_T1w_RAI.1D']);
        filename = fullfile(mridir,'002s01_T1w_RAI.1D'); 
        s01T1wRAI = do_read_coordinates(stimcoor);
        
        figure
        % ft_plot_ortho(funparam, 'transform', data.transform, 'unit', data.unit, 'resolution', resolution, 'style', 'intersect', 'clim', clim);
        ft_plot_ortho(funparam, 'transform', data.transform, 'style', 'intersect');
        axis vis3d
        view([110 36]);
        plot3(fid.pos(NAS_idx,1), fid.pos(NAS_idx,2), fid.pos(NAS_idx,3), 'm.','MarkerSize',80);
        % hold on, plot3(max_x, max_y, max_z, 'm.','MarkerSize',20); % Sub2
        % hold on, plot3(53.6, 1.12, -25.89, 'm.','MarkerSize',20); % Sub2
        hold on, plot3(s01T1wRAI{1}, s01T1wRAI{2}, s01T1wRAI{3}, 'm.','MarkerSize',20);
        
    end
end
