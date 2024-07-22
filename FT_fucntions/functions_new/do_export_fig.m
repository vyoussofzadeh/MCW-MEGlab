function do_export_fig(cfg)

switch cfg.type
    case 'fig'
        filename = [cfg.filename, '.fig'];
        saveas(gcf, fullfile(cfg.outdir, filename));
    case 'pdf'
        filename = [cfg.filename, '.pdf'];
        print(fullfile(cfg.outdirection, filename), '-dpdf');
    case 'png'
        filename = [cfg.filename, '.png'];
        print(fullfile(cfg.outdir, filename), '-dpng', '-r300'); % Save as PNG with 300 DPI resolution
    case 'svg'
        filename = [cfg.filename, '.svg'];
        print(fullfile(cfg.outdir, filename), '-dsvg'); % Save as SVG
end

end
