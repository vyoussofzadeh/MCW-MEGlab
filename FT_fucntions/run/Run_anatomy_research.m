%% Volumetric-based analysis
% mridir = fullfile(indir,subj,'Anatomy/mri');
% d = rdir(fullfile(mridir,subj,'subjectimage*.mat'));

% mridir = fullfile(indir,subj,'Anatomy/mri');
fidfile = fullfile(indir,subj,'Anatomy/bem/Anatomy-fiducials.fif');

% /rcc/stor1/projects/ECP/MEG/MEG_Work/EC1002/Anatomy/mri/T1.mgz
% /rcc/stor1/projects/ECP/MEG/MEG_Work/EC1002/Anatomy/bem/Anatomy-fiducials.fif

clear fid
% if ~isempty(d)
% mripfile = fullfile(mridir,'T1.mgz');
if exist(outputmridir, 'file') == 0, mkdir(outputmridir); end

headshape = ft_read_headshape(datafile);

cfg = [];
cfg.megdata = ep_data.grad;
cfg.mripfile = mripfile;
cfg.hsfile = datafile; % headshape;
cfg.fid = headshape.fid;
cfg.outputmridir = outputmridir;
cfg.subj = subj;
cfg.plotflag = 2;
cfg.atlas = atlas;
cfg.allpath = allpath;
%     [mri_realigned,individual_headmodel,headshape, individual_grid_8mm, individual_grid_10mm] = vy_mri_neuromag2(cfg);
[mri_realigned,individual_headmodel,headshape, individual_grid_8mm, individual_grid_10mm, brain] = vy_mri_neuromag7(cfg);
%     vy_do_freesurfer(cfg);
% end
cd(outd.sub)


%% Choosing mesh
choose_grid = 2;
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