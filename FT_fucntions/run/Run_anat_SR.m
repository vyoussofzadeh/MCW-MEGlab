%% Volumetric-based analysis
bsdir = '/data/MEG/Projects/Cloze2';
mridir = fullfile(bsdir,'anat');
d = rdir(fullfile(mridir,subj,'subjectimage*.mat'));
clear fid
if ~isempty(d)
    sMRI1 = d.name;
    load(sMRI1);
    fid.SCS = SCS;
    fid.NCS = NCS;
    mripfile = fullfile(mridir,subj,'T1.nii');
    if exist(outputmridir, 'file') == 0, mkdir(outputmridir); end
    cfg = [];
    cfg.megdata = cln_data.grad;
    cfg.mripfile = mripfile;
    cfg.hsfile = datafile; % headshape;
    cfg.fid = fid;
    cfg.outputmridir = outputmridir;
    cfg.subj = subj;
    cfg.plotflag = 2;
    cfg.atlas = atlas;
    cfg.indir = indir;
    cfg.outd.sub = outd.sub;
    cfg.flag = flag;
    outanat = vy_mri_neuromag_SR(cfg);
    %     vy_do_freesurfer(cfg);
end
cd(outd.sub)

%% Choosing mesh
choose_grid = 2;
switch choose_grid
    % if flag.warping == 1
    case 1
        switch meshgridres
            case 1
                individual_grid = outanat.individual_grid_10mm_indiv;
            case 2
                individual_grid = outanat.individual_grid_8mm_indiv;
        end
    case 2
        switch meshgridres
            case 1
                meshtag = 'lowres';
                %         load('standard_sourcemodel3d10mm');
                load temp_grid % low-res
                template_grid = ft_convert_units(template_grid, 'mm');
                individual_grid = outanat.individual_grid_10mm;
            case 2
                meshtag = 'highres';
                %         load('standard_sourcemodel3d8mm');
                load temp_grid_8mm % high-res
                individual_grid = outanat.individual_grid_8mm;
        end
end

%% Anatomoy check!
if flag.anatomycheck == 1 
    cfg = [];
    cfg.saveflag = 1;
    cfg.headmodel = outanat.individual_headmodel;
    cfg.leadfield = individual_grid;
    cfg.mri_realigned  = outanat.mri_realigned;
    cfg.headshape = outanat.headshape;
    cfg.outputmridir = outputmridir;
    cfg.mtd = 'vol';
    vy_mri_inspection(cfg, datain);
    %     vy_mri_inspection(t_data, individual_headmodel,individual_grid,headshape, mri_realigned,outputmridir,saveflag);
end

%%
template_mri = ft_read_mri(fullfile(allpath.ft_path,'template/anatomy','single_subj_T1.nii')); %
