%% freq analysis (fft)

% savepath = fullfile(outd.sub,'Freq');
% if exist(savepath, 'file') == 0, mkdir(savepath), end

cfg = [];
cfg.savefile = [];
cfg.saveflag = 0;
cfg.foilim = [2 40];
cfg.plotflag  = 1;
cfg.tapsmofrq       = 4;
cfg.taper    = 'hanning';
do_fft(cfg, datain);

%% Time-freq analysis (tfr)
% savepath = fullfile(outd.sub,'Freq');
% if exist(savepath, 'file') == 0, mkdir(savepath), end

for i=1:length(datain.time)    
    tmax(i) = datain.time{i}(end); 
end

cfg = [];
cfg.savefile = [];
cfg.saveflag = 1;
cfg.lay  = lay;
cfg.subj = subj;
cfg.toi = [datain.time{i}(1),min(tmax)];
tfr = do_tfr(cfg, datain);

%%
% tfr plotting
cfg = [];
cfg.fmax = 40;
cfg.toi = [datain.time{2}(1), min(tmax)];
cfg.savepath = [];
cfg.savefile = 'tfr';
cfg.title = 'test';
cfg.baseline = [-0.3,0];
[time_of_interest,freq_of_interest] = do_tfr_plot(cfg, tfr);

disp(['peaked timing: ',num2str(time_of_interest),' Sec'])
disp(['peaked freq: ',num2str(freq_of_interest),' Hz']);
