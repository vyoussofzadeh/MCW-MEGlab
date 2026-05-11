function [tc, tsec] = pac_aggregate_timecourse(PAC, ip, ia, X, fs, how, K)
% PAC from pac_time_epoch_win; pick band pair (ip, ia)
% how: 'mean'|'median'|'env-weighted'; K is #top channels for weighted (default 20)
Zc = PAC.z{ip,ia};              % C×nbins
tsec = PAC.tbin_sec{ip};        % 1×nb
switch lower(how)
    case 'mean'
        tc = mean(Zc, 1);
    case 'median'
        tc = median(Zc, 1);
    case 'env-weighted'
        if nargin<7 || isempty(K), K=20; end
        env = mean(abs(hilbert(double(X).').'),2);    % C×1
        [~,ix] = maxk(env, min(K, numel(env)));
        w = env(ix) / (sum(env(ix))+eps);
        tc = w.' * Zc(ix,:);                          % 1×nb
    otherwise
        error('how must be mean|median|env-weighted');
end
end
