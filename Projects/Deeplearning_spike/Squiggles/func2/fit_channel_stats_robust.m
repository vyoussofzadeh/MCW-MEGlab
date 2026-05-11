function [med,madv] = fit_channel_stats_robust(C)
% Robust per-channel stats (median, MADsigma) for variable-length epochs.
if isempty(C) || isempty(C{1})
    error('fit_channel_stats_robust:EmptyTraining','XTrain is empty; cannot fit robust stats.');
end
ch = size(C{1},1);
med  = zeros(ch,1);
madv = ones(ch,1);
for r = 1:ch
    rows = cellfun(@(z) double(z(r,:)), C, 'UniformOutput', false);
    rows = rows(~cellfun(@isempty, rows));
    if isempty(rows)
        med(r)  = 0;    % safe fallback
        madv(r) = 1;
        continue;
    end
    x = [rows{:}];                      % 1 x (sum T)
    if isempty(x) || all(isnan(x))
        med(r)  = 0;
        madv(r) = 1;
    else
        m = median(x,'omitnan');
        med(r)  = m;
        madv(r) = 1.4826 * median(abs(x - m),'omitnan');  % MAD?~sigma
        madv(r) = max(madv(r), 1e-12);
    end
end
end
