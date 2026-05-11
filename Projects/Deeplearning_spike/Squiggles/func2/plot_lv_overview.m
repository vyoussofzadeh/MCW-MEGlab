function plot_lv_overview(X, fs)
Xlv = build_longview_features(single(X), fs);
[C,T,~] = size(Xlv);
t = (0:T-1)/fs;

raw = squeeze(Xlv(:,:,1));
gfp = std(raw,[],1);

figure('Name','LV overview','Color','w');

subplot(4,1,1)
plot(t, gfp, 'k','LineWidth',1.2); grid on
title('GFP (raw)'); xlim([t(1) t(end)])

subplot(4,1,2)
imagesc(t,1:C, squeeze(Xlv(:,:,3))); axis tight xy; colormap(parula); colorbar
title('amp (relative)'); ylabel('Channel')

subplot(4,1,3)
imagesc(t,1:C, squeeze(Xlv(:,:,5))); axis tight xy; colormap(parula); colorbar
title('slope'); ylabel('Channel')

subplot(4,1,4)
imagesc(t,1:C, squeeze(Xlv(:,:,7))); axis tight xy; colormap(parula); colorbar
title('sharpness'); ylabel('Channel'); xlabel('Time (s)')
end
