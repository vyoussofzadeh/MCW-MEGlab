function plot_feature_maps_epoch(Xstack, fs, names)
% Xstack: C×T×D
if nargin<3 || isempty(names), names = {'raw','d1','d2','env','smooth','edge'}; end
D = size(Xstack,3); t = (0:size(Xstack,2)-1)/fs;
figure('Color','w'); tiledlayout(D,1,'TileSpacing','compact','Padding','compact');
for k=1:D
    nexttile; imagesc(t,1:size(Xstack,1), Xstack(:,:,k));
    axis tight xy; colormap(parula); colorbar; title(names{k});
    if k<D, set(gca,'XTickLabel',[]); else, xlabel('Time (s)'); ylabel('Chan'); end
end
end

% Example:
% plot_feature_maps_epoch(Xtr2D{1}, fs, {'raw','d1','d2','env'});
