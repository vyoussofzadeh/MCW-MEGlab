function do_export_ft_sourcemaps(cfg)

switch cfg.type
    case 'fig'
        filename = [cfg.filename,'.fig'];
        saveas(gcf, fullfile(cfg.outdir,filename));
    case 'pdf'
        filename = [cfg.filename,'.pdf'];
        print(fullfile(cfg.outdir,filename), '-dpdf');
end