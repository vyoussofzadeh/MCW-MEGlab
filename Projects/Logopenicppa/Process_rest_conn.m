%% Logopenicppa dataset, Medical College of Wisconsin

% Script: BS Process (seed-based conn analysis)
% Project: Logopenicppa_rest
% Writtern by: Vahab YoussofZadeh
% Update: 04/12/2023

clear; clc, close('all'); warning off,

%% Paths
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/Logopenicppa/functions')
[atlas, all_path] = do_setpath();

%%
indir = '/data/MEG/Research/logopenicppa';
outdir = fullfile(indir,'ft_process');
rawdatadir = fullfile(indir,'raw');

%% Analyses
flag = [];
flag.preproces.filtering = 1;
flag.preproces.artifact = 1;
flag.preproces.ica = 1;
flag.ic_selection = 1;
flag.seed_check = 1;
flag.anatomy_check = 2;

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
datafile.subj

%%
cfg = []; cfg.layout = 'neuromag306mag.lay'; lay = ft_prepare_layout(cfg);

%%
mridir = '/data/MEG/Research/logopenicppa/MRI/idsc9005/stimulation_rois';

%%
for i=1:length(datafile.task_run)
    close all
    datafile_sel = datafile.datafile_fif{i}; % spm_select(inf,'dir','Select MEG folder'); % e.g. H:\VNS\MEG\C-105\CRM\1
    Index = strfind(datafile_sel, '/');
    subj = datafile_sel(Index(6)+1:Index(7)-3);
    Index = strfind(datafile_sel, '32037');
    Index1 = strfind(datafile_sel, 'Run');
    run = datafile_sel(Index1+4);
    %     if isempty(Index)
    %         Index = strfind(datafile_sel, 'Run');
    %     end
    %     if isempty(Index)
    %         break,
    %     else
    session  = datafile_sel(Index+10);
    tkz = tokenize(subj,'_');
    mripfile = fullfile(mridir,[tkz{2}, 's01_T1w.nii']);
    subdir = fullfile(outdir,subj, tag, session);
    if exist(subdir, 'file') == 0
        mkdir(subdir);   %create the directory
    end
    cd(subdir)
    disp(['outputdir:',subdir])
    
%     savefile_seed_conn = fullfile(subdir,['seed_conn_',subj,'_session_', session, '_Run', run, '.mat']);
    savefile_seed_conn = fullfile(subdir,['seed_conn_',subj,'_run_', run, '.mat']); 
    if exist(mripfile,'file')== 2 && exist(savefile_seed_conn,'file') ~=2
        
        disp(datafile_sel)
        disp(['subj:',subj])
        disp(['Run:',run])
        
        %-elec/grad
        sens = ft_read_sens(datafile_sel);
        sens = ft_convert_units(sens,'mm');

        %% Updating datalog
        Datalog = [];
        Datalog.subj = subj;
        Datalog.run = run;
        Datalog.datafile = datafile;
        Datalog.subdir = subdir;
        Datalog.outdir = outdir;
        Datalog.mripfile = mripfile;
        Datalog.mridir = mridir;
        
        %% Preprocesssing
        cfg = [];
        cfg.datalog = Datalog;
        cfg.flag = flag;
        cfg.lay = lay;
        cfg.all_path = all_path;
        cln_data = do_preprocess_rest(cfg, datafile_sel);
        
        %% setpath again
        addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/Logopenicppa/functions')
        [atlas, all_path] = do_setpath();
        
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
        cfg.all_path = all_path;
        anat = do_anatomy(cfg);
        
        %%
        individual_headmodel = anat.individual_headmodel;
        mri_realigned = anat.mri_realigned;
        headshape = anat.headshape;
        
        %% Choosing mesh
        flag.meshgrid_sel = 1;
        switch flag.meshgrid_sel
            case 1
                meshtag = 'lowres';
                load temp_grid % low-res
                template_grid = ft_convert_units(template_grid, 'mm');
                individual_grid = anat.individual_grid_10mm;
            case 2
                meshtag = 'highres';
                load temp_grid_8mm % high-res
                individual_grid = anat.individual_grid_8mm;
        end
        
        %% Anatomoy check!
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
        anat.individual_grid = individual_grid;
        anat.template_grid = template_grid;
        anat.atlas = atlas;
        anat.sens = sens;
        
        %%
        cfg = [];
        cfg.anat = anat;
        cfg.Datalog = Datalog;
        seed_coor = do_seed_inspection(cfg);
        
        %%    
        pflag = [];
        pflag.allconn = 2;
        pflag.grid = 2;
        pflag.grid_seed = 2;
        pflag.aal = 2;
        
        sflag = [];
        sflag.seed = 1;
        sflag.seed_map = 1;
        sflag.seed_conn = 1;
 
        cd(subdir)
        cfg = [];
        cfg.cov_matrix = cov_matrix;
        cfg.seed = seed_coor;
        cfg.anat = anat;
        cfg.foi = [18,25];
        cfg.pflag = pflag;
        cfg.sflag = sflag;
        net_conn_seed = do_conn_seed(cfg, cln_data);
        
        %%
%         cfg = [];
%         cfg.cov_matrix = cov_matrix;
%         cfg.anat = anat;
%         cfg.foi = net_conn_seed.foi;
%         net_conn = do_wPLIconn1(cfg, net_conn_seed.source_active);
        
        %%
        save(savefile_seed_conn,'net_conn_seed');
%         savefile_whole_conn = fullfile(subdir,['wholeb_conn_',subj,'_run_', run, '.mat']); save(savefile_whole_conn,'net_conn')
    end
end
