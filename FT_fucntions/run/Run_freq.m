%% freq analysis (fft)

savepath = fullfile(outd.sub,'Freq');
if exist(savepath, 'file') == 0, mkdir(savepath), end

cfg = [];
cfg.savefile = fullfile(savepath,[dtag.dcon,'_fft_',subj,'.mat']);
cfg.saveflag = 1;
cfg.foilim = [2 40];
cfg.plotflag  = 1;
cfg.tapsmofrq       = 4;
cfg.taper    = 'hanning';
vy_fft(cfg, datain);

%% Time-freq analysis (tfr)
savepath = fullfile(outd.sub,'Freq');
if exist(savepath, 'file') == 0, mkdir(savepath), end

tmax = [];
for i=1:length(datain.time)    
    tmax(i) = datain.time{i}(end); 
end

cfg = [];
cfg.savefile = fullfile(savepath,[dtag.dcon,'_tfr_',subj,'.mat']);
cfg.saveflag = 1;
cfg.lay  = lay;
cfg.subj = subj;
cfg.toi = [datain.time{i}(1),min(tmax)];
tfr = vy_tfr(cfg, datain);

%%
% tfr plotting
cfg = [];
cfg.baseline = [-1 0];
cfg.fmax = fmax;
cfg.toi = [datain.time{2}(1), min(tmax)];
cfg.savepath = fullfile(outd.sub,'Freq');
cfg.savefile = fullfile(savepath,[dtag.dcon,'_tfr2_',subj]);
[time_of_interest,freq_of_interest] = vy_tfr_plot(cfg, tfr);

disp(['peaked timing: ',num2str(time_of_interest),' Sec'])
disp(['peaked freq: ',num2str(freq_of_interest),' Hz']);
