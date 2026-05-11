function plot_pac_tc_avg(X, PAC, ip, ia, fs, K)
if nargin<6, K=20; end
env = abs(hilbert(double(X).').');        % C×T
mEnv = mean(env,2);
[~,ix] = maxk(mEnv, min(K, size(X,1)));
Z = PAC.z{ip,ia};  % C×nb
t = PAC.tbins{ip}/fs;
figure('Color','w'); plot(t, mean(Z(ix,:),1), 'LineWidth',1.4); grid on
xlabel('Time (s)'); ylabel('zPAC'); title(sprintf('Mean zPAC over top-%d sensors',numel(ix)));
end
