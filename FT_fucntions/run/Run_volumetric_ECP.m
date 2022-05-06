%% Volumetric-based analysis
% mridir = fullfile(indir,subj,'Anatomy/mri');
% d = rdir(fullfile(mridir,subj,'subjectimage*.mat'));

mridir = fullfile(indir,subj,'Anatomy/mri');
fidfile = fullfile(indir,subj,'Anatomy/bem/Anatomy-fiducials.fif');

% /rcc/stor1/projects/ECP/MEG/MEG_Work/EC1002/Anatomy/mri/T1.mgz
% /rcc/stor1/projects/ECP/MEG/MEG_Work/EC1002/Anatomy/bem/Anatomy-fiducials.fif

clear fid
% if ~isempty(d)
mripfile = fullfile(mridir,'T1.mgz');
if exist(outputmridir, 'file') == 0, mkdir(outputmridir); end

headshape = ft_read_headshape(datafile);

cfg = [];
cfg.megdata = ep_data.app.grad;
cfg.mripfile = mripfile;
cfg.hsfile = datafile; % headshape;
cfg.fid = headshape.fid;
cfg.outputmridir = outputmridir;
cfg.subj = subj;
cfg.plotflag = 2;
cfg.atlas = atlas;
%     [mri_realigned,individual_headmodel,headshape, individual_grid_8mm, individual_grid_10mm] = vy_mri_neuromag2(cfg);
[mri_realigned, individual_headmodel, headshape, individual_grid_8mm, individual_grid_10mm] = vy_mri_neuromag5(cfg);
%     vy_do_freesurfer(cfg);
% end
cd(outd.sub)


%% Choosing mesh
switch choose_grid
    % if flag.warping == 1
    case 1
        switch flag.meshgrid_sel
            case 1
                individual_grid = outanat.individual_grid_10mm_indiv;
            case 2
                individual_grid = outanat.individual_grid_8mm_indiv;
        end
    case 2
        switch flag.meshgrid_sel
            case 1
                meshtag = 'lowres';
                %         load('standard_sourcemodel3d10mm');
                load temp_grid % low-res
                template_grid = ft_convert_units(template_grid, 'mm');
                individual_grid = individual_grid_10mm;
            case 2
                meshtag = 'highres';
                %         load('standard_sourcemodel3d8mm');
                load temp_grid_8mm % high-res
                individual_grid = individual_grid_8mm;
        end
        % else
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
    vy_mri_inspection(cfg, ep_data);
    %     vy_mri_inspection(t_data, individual_headmodel,individual_grid,headshape, mri_realigned,outputmridir,saveflag);
end

%%
template_mri = ft_read_mri(fullfile(allpath.ft_path,'template/anatomy','single_subj_T1.nii')); %

%%
switch method
    case 1
        %% LCMV
        %         mtag = 'lcmv_stat';
        mtag = 'lcmv';
        outd.vol = fullfile(outd.sub, mtag);
        cfg = [];
        cfg.allpath = allpath;
        cfg.grid = individual_grid;
        cfg.headmodel = individual_headmodel;
        cfg.subj = subj;
        cfg.sens = sens;
        %         cfg.mtag = 'lcmv'; cfg.filterflag =  1;
        cfg.mtag = mtag; cfg.filterflag =  1;
        outd.vol = fullfile(outd.sub,cfg.mtag);
        %         cfg.fb   = [12,30]; % Hz
        %         cfg.fb   = [16,25]; % Hz
        cfg.fb   = [14,27]; % Hz
        cfg.atlas = allpath.atlas_path;
        cfg.outputdir = outd.vol;
        cfg.toi       = {toi};
        cfg.template_grid = template_grid;
        cfg.template_mri = template_mri;
        cfg.plotflag     = 1;
        switch mtag
            case 'lcmv'
                vy_source_lcmv(cfg, cln_data);
            case 'lcmv_stat'
                vy_source_lcmv_stats(cfg, cln_data);
        end
        
    case 2
        %%
        cfg = [];
        mtag = 'conn';
        %         mtag = 'conn_bl'; cfg.fl = [1, 10];% band limited
        %         mtag = 'conn_bs'; cfg.fb = 10; % band-stop
        outd.vol = fullfile(outd.sub,mtag);
        cfg.allpath = allpath;
        cfg.grid = individual_grid;
        cfg.headmodel = individual_headmodel;
        cfg.subj = subj;
        cfg.sens = sens;
        cfg.mtag = mtag;
        cfg.toi  = toi;
        cfg.atlas = allpath.atlas_path;
        cfg.outputdir = outd.vol;
        cfg.template_grid = template_grid;
        cfg.template_mri = template_mri;
        switch mtag
            case 'conn'
                vy_network_light1(cfg,ep_data) % conn-network analysis
            case {'conn_bs','conn_bl'}
                vy_network_light1(cfg, cln_data) % conn-network analysis
        end
        
    case 3
        %%
        disp('1: DICS, simple contrasting');
        disp('2: DICS-stats');
        %          mtag = input(':');
        mtag = 1;
        
        switch mtag
            case 1
                mtag = 'dics';
            case 2
                mtag = 'dics_stat';
        end
        
        %         mtag = 'dics_ratio';
        %         mtag = 'dics_fs';
        if rl == 1, mtag_lab = 'dics_res'; else, mtag_lab = 'dics'; end
        outd.vol = fullfile(outd.sub,mtag_lab);
        
        switch mtag
            case 'dics'
                cfg = [];
                cfg.grid = individual_grid;
                cfg.allpath = allpath;
                cfg.freq_of_interest  = freq_of_interest; % Hz
                cfg.headmodel = individual_headmodel;
                cfg.sens = sens;
                cfg.mtag = mtag;
                cfg.subj = subj;
                cfg.toi = toi;
                cfg.outputdir = outd.vol;
                switch choose_grid
                    case 1
                        cfg.template_grid = [];
                        cfg.mri_aligned = outanat.mri_realigned;
                    case 2
                        cfg.template_grid = template_grid;
                end
                cfg.template_mri = template_mri;
                cfg.fmax = fmax;
                cfg.savedata = fullfile(outd.vol,[mtag,'_',subj]);
                cfg.flag = flag;
                cfg.plotting = 1;              
                vy_source_dics(cfg, ep_data);
                
            case 'dics_ratio'
                
                mtag = 'dics_ratio';
                cfg = [];
                cfg.grid = individual_grid;
                cfg.allpath = allpath;
                cfg.freq_of_interest   = freq_of_interest; % Hz
                cfg.headmodel = individual_headmodel;
                cfg.sens = sens;
                cfg.mtag = mtag;
                cfg.subj = subj;
                cfg.toi = toi;
                cfg.outputdir = outd.vol;
                cfg.template_grid = template_grid;
                cfg.template_mri = template_mri;
                cfg.savedata = fullfile(outd.vol,[mtag,'_',subj,'.mat']);
                vy_source_dics_ratio(cfg, ep_data);
                
            case 'dics_fs'
                %%
                anatomy_dir     = '/data/MEG/Clinical/ft_process/19/bednar_peggy/anat/BAK';
                load(fullfile(anatomy_dir,[subj,'_headmodel.mat']));
                load(fullfile(anatomy_dir,[subj,'_leadfield.mat']));
                load(fullfile(anatomy_dir,[subj,'_sourcemodel.mat']));
                
                if anatomy_check_flag == 1
                    cfg1 = [];
                    cfg1.saveflag = 2;
                    cfg1.headmodel = headmodel;
                    cfg1.sourcemodel = sourcemodel;
                    cfg1.leadfield = leadfield;
                    cfg1.mri = anatomy_dir;
                    cfg1.mtd = 'surf';
                    cfg1.headshape = headshape;
                    vy_mri_inspection(cfg1, t_data);
                    %                 individual_headmodel,individual_grid,headshape, mri_realigned,outputmridir,saveflag
                end
                cfg = [];
                cfg.headmodel = headmodel;
                cfg.sourcemodel = sourcemodel;
                cfg.leadfield   = leadfield;
                cfg.mtag = mtag;
                cfg.sens = sens;
                cfg.subj = subj;
                cfg.outputdir = outd.vol;
                vy_source_dics_fs(cfg, ep_data);
                
            case 'dics_stat'
                mtag = 'dics_stat';
                cfg = [];
                cfg.grid = individual_grid;
                cfg.allpath = allpath;
                %                 cfg.f   = 20; % Hz
                cfg.headmodel = individual_headmodel;
                cfg.sens = sens;
                cfg.mtag = mtag;
                cfg.subj = subj;
                cfg.toi = toi;
                cfg.outputdir = outd.vol;
                cfg.template_grid = template_grid;
                cfg.template_mri = template_mri;
                cfg.flag = flag;
                cfg.savedata = fullfile(outd.vol,[mtag,'_',subj,'.mat']);
                vy_source_dics_stats(cfg, ep_data);
        end
        %% mne_ica
    case 4
        switch dc
            case 1
                datain = cln_data_str;
            case 2
                datain = cln_data_math;
        end
        Run_icaMNE
        
end
