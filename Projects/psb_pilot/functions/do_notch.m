function dataout = do_notch(cfg_main, datain)


cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.foilim = [2 100];
cfg.plotflag  = 1;
cfg.tapsmofrq = 8;
cfg.taper     = 'hanning';
[freq,ff,psd] = do_fft(cfg, datain);
grid on
grid minor
title(['Before band-stop filtering-',cfg_main.subj]);

%- finding notch freq
idx = find(freq.freq ==30); TF = islocalmax(psd); TF(1:idx-5) = 0;
hold on
plot(ff(TF),psd(TF),'r*')

idx2 = find(TF == 1);
if length(idx2)==1
    fsb = round(ff(idx2));
else
    [val, idx] = max(psd(TF));
    fff = ff(TF);
    fsb = round(fff(idx));
end
disp([num2str(fsb),'Hz freq was selected for notch filteting']);

%% updated, 02/24/22
cfg = [];
cfg.bsfilter = 'yes';
%     cfg.bsfreq = [29 32]; % or whatever you deem appropriate
cfg.bsfreq = [fsb-1 fsb+1]; % or whatever you deem appropriate
dataout = ft_preprocessing(cfg, datain);

%%
cfg = [];
cfg.savefile = [];
cfg.saveflag = 2;
cfg.foilim = [2 100];
cfg.plotflag  = 1;
cfg.tapsmofrq = 5;
cfg.taper     = 'hanning';
do_fft(cfg, dataout);
grid on
grid minor
title(['After band-stop filtering-',cfg_main.subj]);

end