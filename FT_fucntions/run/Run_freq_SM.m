%% freq analysis (fft)

% savepath = fullfile(outd.sub,'Freq');
% if exist(savepath, 'file') == 0, mkdir(savepath), end

% cfg = [];
% cfg.savefile = fullfile(savepath,[stag,'_fft_',subj,'.mat']);
% cfg.saveflag = 1;
% cfg.foilim = [2 40];
% cfg.plotflag  = 1;
% cfg.tapsmofrq       = 4;
% cfg.taper    = 'hanning';
% vy_fft(cfg, datain);

%% Time-freq analysis (tfr)
savepath = fullfile(outd.sub,'Freq');
if exist(savepath, 'file') == 0, mkdir(savepath), end
cfg = [];
cfg.savefile = savetag;
cfg.saveflag = 1;
cfg.lay  = lay;
cfg.subj = subj;
tfr = vy_tfr(cfg, datain);

%%
% tfr plotting
cfg = [];
cfg.savepath = fullfile(outd.sub,'Freq');
cfg.savefile = fullfile(savepath,[stag,'_tfr2_',subj]);
[time_of_interest,freq_of_interest] = vy_tfr_plot(cfg, tfr);

disp(['peaked timing: ',num2str(time_of_interest),' Sec'])
disp(['peaked freq: ',num2str(freq_of_interest),' Hz']);
