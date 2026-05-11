function plot_lv_channel(X, fs, chIdx)
% X: [C x T] (epoch), fs in Hz, chIdx (optional)
Xlv = build_longview_features(single(X), fs);     % [C x T x 7]
[C,T,~] = size(Xlv);
t = (0:T-1)/fs;

% pick channel
if nargin<3 || isempty(chIdx)
    [~, chIdx] = max(rms(X,2)); % most energetic sensor
end

names = {'raw','topo','amp','meanAmp','slope','halfSlope','sharpness'};
figure('Name',sprintf('LV features  ch %d', chIdx), 'Color','w');
for k = 1:7
    subplot(7,1,k)
    plot(t, squeeze(Xlv(chIdx,:,k)), 'LineWidth', 1.0)
    xlim([t(1) t(end)]); grid on
    ylabel(names{k})
    if k==1, title(sprintf('Channel %d (fs=%g Hz)', chIdx, fs)); end
end
xlabel('Time (s)')
end
