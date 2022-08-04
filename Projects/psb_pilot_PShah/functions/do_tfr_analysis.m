function [toi,foi, tfr]  = do_tfr_analysis(cfg_main, datain)

tmax = []; 
for i=1:length(datain.time)
    tmax(i) = datain.time{i}(end); 
end

cfg = [];
cfg.savefile = []; 
cfg.saveflag = 1; 
cfg.lay  = cfg_main.layout;
cfg.subj = num2str(cfg_main.subj); 
cfg.toi = [datain.time{i}(1),min(tmax)];
tfr = do_tfr(cfg, datain);

cfg = []; 
cfg.baseline = [cfg_main.baseline(1) cfg_main.baseline(2)]; 
cfg.fmax = cfg_main.fmax; 
cfg.toi = [datain.time{2}(1), min(tmax)]; 
cfg.savepath = []; 
cfg.savefile = [];
cfg.title = cfg_main.title;
[toi,foi] = do_tfr_plot(cfg, tfr);

end