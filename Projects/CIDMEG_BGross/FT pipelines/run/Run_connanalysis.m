% pause, close all,
cfg            = [];
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.foilim     = [1 40];
cfg.tapsmofrq  = 2;
cfg.keeptrials = 'yes';
cfg.pad = 4;
freq    = ft_freqanalysis(cfg, vs_roi1);

cfg         = [];
cfg.method    = 'wpli_debiased';
source_conn = ft_connectivityanalysis(cfg, freq);
par = 'wpli_debiasedspctrm';

% cfg         = [];
% cfg.method    = 'amplcorr';
% source_conn = ft_connectivityanalysis(cfg, freq);
% par = 'amplcorrspctrm';

%%
source_conn1 = source_conn;
for i=1:size(source_conn1.(par),3)
    tmp = source_conn1.(par);
    tmp(:,:,i) = tril(squeeze(tmp(:,:,i)));
end
source_conn1.(par) = tmp;

n = floor(size(source_conn1.(par),1)/15);

y = round(linspace(1,size(source_conn1.(par),1),n));

figure
cfg           = [];
cfg.parameter = par;
cfg.zlim      = [0 1];
cfg.channel = y(1):y(1+1);
ft_connectivityplot(cfg, source_conn1);
set(gcf, 'Position', [800   400   900   900]);
%             pause,

