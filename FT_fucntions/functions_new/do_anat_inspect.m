function do_anat_inspect(cfg, ~)


switch cfg.mtd
    case 'vol'
        %%
        figure;
        ft_plot_vol(cfg.headmodel, 'facecolor', 'cortex', 'edgecolor', 'none');alpha 0.5; camlight;
        hold on;
        ft_plot_headshape(cfg.headshape);
        ft_plot_mesh(cfg.leadfield.pos(cfg.leadfield.inside, :));
        view ([0 90])
        if ~isempty(cfg.saveflag)
            savepath = fullfile(cfg.outputmridir,'headshape');
            hcp_write_figure([savepath,'.png'], gcf, 'resolution', 300);
        end
        
    case 'surf'
        %%
        figure; hold on;
        ft_plot_vol(cfg.headmodel, 'facecolor', 'none'); alpha 0.5;
        ft_plot_mesh(cfg.sourcemodel, 'facecolor', 'cortex', 'edgecolor', 'none'); camlight;
        
        %%
        figure; hold on;
        ft_plot_mesh(cfg.sourcemodel, 'facecolor', 'cortex', 'edgecolor', 'none'); camlight;
        ft_plot_mesh(cfg.leadfield.pos(cfg.leadfield.inside, :));
        
end