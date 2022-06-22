antdir = fullfile(outdir, subj,'anat'); % output dir
if exist(antdir, 'file') == 0, mkdir(antdir); end
cd(antdir)

mripfile = fullfile(mridir,mrifile);
headshape = ft_read_headshape(datafile);

cfg = [];
cfg.megdata = data_clean.grad;
cfg.mripfile = mripfile;
cfg.hsfile = datafile; % headshape;
cfg.fid = headshape.fid;
cfg.outputmridir = antdir;
cfg.subj = subj;
cfg.plotflag = 1;
[mri_realigned,individual_headmodel,headshape, individual_grid_8mm, individual_grid_10mm] = do_anatomy(cfg);

%- Choosing mesh (as needed for source modelling)
flag.meshgrid_sel = 1;
choose_grid = 2;
switch choose_grid
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