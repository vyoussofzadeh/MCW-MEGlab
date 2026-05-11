function plot_lv_heatmaps(X, fs)
Xlv = build_longview_features(single(X), fs);   % [C x T x 7]
[C,T,~] = size(Xlv);
t = (0:T-1)/fs;
names = {'raw','topo','amp','meanAmp','slope','halfSlope','sharpness'};

figure('Name','LV maps (channel × time)', 'Color','w');
for k = 1:7
    subplot(3,3,k)
    imagesc(t, 1:C, squeeze(Xlv(:,:,k)))
    axis tight xy
    colormap(parula); colorbar
    title(names{k})
    xlabel('Time (s)'); ylabel('Channel')
end
end
