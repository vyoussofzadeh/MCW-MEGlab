%% Volumetric-based analysis
% mridir = fullfile(indir,subj,'brainstorm_db/anat');
mridir = '/MEG_data/Vahab/Projects/Oddball/anat';
% d = rdir(fullfile(mridir,subj,'subjectimage*.mat'));
mripfile = fullfile(mridir, 'sub-RV1217_ses-01_T1w_1.nii');

%%
individual_mri = ft_read_mri(mripfile);
ft_sourceplot([], individual_mri);
individual_mri = ft_convert_units(individual_mri, 'mm');

%% Normalized coordinate system (NCS)
% cfg          = [];
% % cfg.method   = 'interactive';
% cfg.method = 'fiducial'; % the following voxel coords were determined interactive
% cfg.coordsys = 'ctf';
% cfg.fiducial.ac   = fid.NCS.AC;
% cfg.fiducial.pc   = fid.NCS.PC;
% cfg.fiducial.xzpoint  = fid.NCS.IH;
% cfg.fiducial.right   = fid.SCS.RPA;
% mri_realigned     = ft_volumerealign(cfg, individual_mri);

% nas = [-8.673617379884035e-19, 0.0877611252118013, 0.0];
% rpa = [0.06523973384436928, -3.8030137245920384e-20, 0.0];
% lpa = [-0.0725013953618283, -4.716869082909292e-19, 0.0];

nas = [-7.949005, -115.6568, -1.4552];
lpa = [66.451, -37.60901, -37.315];
rpa = [-71.54901, -33.39021, -32.0415];


rpa                = [-nas(2) nas(1) nas(3)];
nas                = [-lpa(2) lpa(1) lpa(3)];
lpa                = [-rpa(2) rpa(1) rpa(3)];


% nas(1) = -nas(1); 


% headshape = ft_read_headshape(datafile);
% nas = headshape.fid.pos(1,:);
% lpa = headshape.fid.pos(2,:);
% rpa = headshape.fid.pos(3,:);

% nas =   [  8.8600         0         0];
% lpa =  [  -0.2635    6.8660         0];
% rpa = [    0.2635   -6.8660         0];
% 
%%
% cfg  = [];
% individual_mri = ft_volumereslice(cfg,individual_mri);

%% Subject Coordinate System (SCS / CTF)
cfg = [];
cfg.method = 'fiducial';
cfg.coordsys = 'ctf';
cfg.fiducial.nas   = nas;
cfg.fiducial.lpa   = lpa;
cfg.fiducial.rpa   = rpa;
% cfg.method   = 'interactive';
cfg.spmversion     = 'spm12';
mri_realigned = ft_volumerealign(cfg, individual_mri);

%%
% ft_sourceplot([], mri_realigned);


%%
ft_sourceplot([], mri_realigned); hold on
% figure, hold on
plot3(nas(1,1), nas(1,2), nas(1,3), 'm*');
plot3(lpa(1,1), lpa(1,2), lpa(1,3), 'm*');
plot3(rpa(1,1), rpa(1,2), rpa(1,3), 'm*');
% %%

%%
individual_mri.coordsys = 'spm';
cfg = [];
cfg.output = 'brain';
cfg.spmversion = 'spm12';
individual_seg = ft_volumesegment(cfg, individual_mri);

%%
cfg = [];
cfg.method = 'singleshell';
cfg.spmversion = 'spm12';
individual_headmodel = ft_prepare_headmodel(cfg, individual_seg);

%%
cfg = [];
cfg.method = 'singleshell';
cfg.spmversion = 'spm12';
individual_headmodel = ft_prepare_headmodel(cfg, mri);



%% Source model, warpping with template
load temp_grid % low-res
% load('standard_sourcemodel3d10mm');sourcemodel = ft_convert_units(sourcemodel, 'mm');
cfg                 = [];
cfg.warpmni    = 'yes';
cfg.spmversion     = 'SPM12';
cfg.grid.nonlinear  = 'yes';
cfg.grid.template   = template_grid;
% cfg.grid.template   = sourcemodel;
cfg.mri             = mri;
cfg.grid.unit       = 'mm';
individual_grid_10mm     = ft_prepare_sourcemodel(cfg);

%%
vol = ft_convert_units(vol, 'mm');

figure;
ft_plot_vol(vol, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
hold on;
% ft_plot_headshape(headshape);
ft_plot_mesh(individual_grid_10mm.pos(individual_grid_10mm.inside, :));
view ([0 90])
    

%%
% clear fid
% if ~isempty(d)
%     sMRI1 = d.name;
%     load(sMRI1);
%     fid.SCS = SCS;
%     fid.NCS = NCS;
%     mripfile = fullfile(mridir,'T1.nii');
%     if exist(outputmridir, 'file') == 0, mkdir(outputmridir); end
%     cfg = [];
%     cfg.megdata = t_data.pst.grad;
%     cfg.mripfile = mripfile;
%     cfg.hsfile = datafile; % headshape;
%     cfg.fid = fid;
%     cfg.outputmridir = outputmridir;
%     cfg.subj = subj;
%     cfg.plotflag = 2;
%     cfg.atlas = atlas;
%     cfg.indir = indir;
%     cfg.outd.sub = outd.sub;
%     cfg.flag = flag;
%     outanat = vy_mri_neuromag2(cfg);
%     %     vy_do_freesurfer(cfg);
% end

%% Choosing mesh
% if flag.warping == 1
%     switch meshgrid
%         case 1
%             meshtag = 'lowres';
%             %         load('standard_sourcemodel3d10mm');
%             load temp_grid % low-res
%             template_grid = ft_convert_units(template_grid, 'mm');
%             individual_grid = outanat.individual_grid_10mm;
%         case 2
%             meshtag = 'highres';
%             %         load('standard_sourcemodel3d8mm');
%             load temp_grid_8mm % high-res
%             individual_grid = outanat.individual_grid_8mm;
%     end
% else
%     switch meshgrid
%         case 1            
%             individual_grid = outanat.individual_grid_10mm;
%         case 2
%             individual_grid = outanat.individual_grid_8mm;
%     end
% end

%% Anatomoy check!
saveflag = 2;
if flag.anatomycheck == 1
    
    cfg = [];
    cfg.saveflag = 1;
    cfg.headmodel = vol;
    cfg.leadfield = individual_grid_10mm;
    cfg.mri_realigned  = mri;
%     cfg.headshape = outanat.headshape;
    cfg.outputmridir = './anat';
    cfg.mtd = 'vol';
    vy_mri_inspection(cfg, t_data);
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
        cfg.template_mri = m;
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
        cfg.headmodel = outanat.individual_headmodel;
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
                vy_network_light1(cfg,t_data) % conn-network analysis
            case {'conn_bs','conn_bl'}
                vy_network_light1(cfg, cln_data) % conn-network analysis
        end
        
    case 3
        %%
        mtag = 'dics';
%         mtag = 'dics_ratio';
%                 mtag = 'dics_stat';
        %         mtag = 'dics_fs';
        outd.vol = fullfile(outd.sub,mtag);
        
        switch mtag
            case 'dics'
                cfg = [];
                cfg.grid = individual_grid_10mm;
                cfg.allpath = allpath;
                cfg.freq_of_interest  = freq_of_interest; % Hz
                cfg.headmodel = vol;
                cfg.sens = sens;
                cfg.mtag = mtag;
                cfg.subj = subj;
                cfg.toi = toi;
                cfg.outputdir = './';
                if flag.warping ==1
                cfg.template_grid = template_grid;
                end
                cfg.template_mri = mri;
                cfg.fmax = fmax;
                cfg.savedata = './fq';
                cfg.flag = flag;
                vy_source_dics(cfg, ep_data);
            case 'dics_ratio'
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
                cfg = [];
                cfg.grid = individual_grid;
                cfg.allpath = allpath;
                cfg.f   = 20; % Hz
                cfg.headmodel = individual_headmodel;
                cfg.sens = sens;
                cfg.mtag = mtag;
                cfg.subj = subj;
                cfg.toi = toi;
                cfg.outputdir = outd.vol;
                cfg.template_grid = template_grid;
                cfg.template_mri = template_mri;
                vy_source_dics_stats(cfg, ep_data);
        end
        
end
