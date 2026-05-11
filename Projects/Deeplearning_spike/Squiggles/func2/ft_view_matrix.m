% X: [channels x time], fs: sampling rate (Hz), labels: cellstr of channel names
function ft_view_matrix(X, fs, labels)
if nargin<3 || isempty(labels)
    labels = arrayfun(@(i) sprintf('Ch%03d',i), 1:size(X,1), 'uni', 0);
end

data = [];
data.fsample = fs;
data.label   = labels(:);
data.trial   = { double(X) };                    % FieldTrip expects double
data.time    = { (0:size(X,2)-1)/fs };

cfg = [];
cfg.viewmode  = 'vertical';                      % one channel per row
cfg.blocksize = min(5, size(X,2)/fs);            % seconds per page
cfg.ylim      = 'maxmin';                        % autoscale per channel
cfg.preproc.demean = 'yes';                      % quick baseline
% optional:
% cfg.channel  = {'MEG*'};                       % select subset
% cfg.zoomsig  = 1.5;                            % vertical zoom
% cfg.artfctdef.clip = 'yes';                    % mark clips

ft_databrowser(cfg, data);                       % (aka ft_databrowse in some versions)
end