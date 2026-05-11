function plot_epoch_maps(X2D_epoch, fs, names)
if nargin<3, names=[]; end
D = size(X2D_epoch,3); t=(0:size(X2D_epoch,2)-1)/fs;
figure('Color','w'); tiledlayout(D,1,'TileSpacing','compact','Padding','compact');
for k=1:D
    nexttile; imagesc(t,1:size(X2D_epoch,1), X2D_epoch(:,:,k));
    axis tight xy; colormap(parula); colorbar
    if ~isempty(names) && k<=numel(names), title(names{k}); else, title(sprintf('map %d',k)); end
    if k<D, set(gca,'XTickLabel',[]); else, xlabel('Time (s)'); ylabel('Channel'); end
end
end
