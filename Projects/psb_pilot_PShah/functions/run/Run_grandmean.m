a_data = do_ave(datain);
% savepath = fullfile(outd.sub,'Timelock');
% if exist(savepath, 'file') == 0, mkdir(savepath), end
cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.lay  = lay;
do_ave_plot(cfg, a_data);