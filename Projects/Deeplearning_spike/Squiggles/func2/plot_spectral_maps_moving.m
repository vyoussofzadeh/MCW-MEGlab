function plot_spectral_maps_moving(B, fs)
% B from feats_spectral_moving: fields bp1..bpK (C×T)
K = numel(fieldnames(B)); fn = fieldnames(B);
t = (0:size(B.bp1,2)-1)/fs; C = size(B.bp1,1);
figure('Color','w'); tiledlayout(K,1,'TileSpacing','compact','Padding','compact');
for k=1:K
    nexttile; M = B.(fn{k}); imagesc(t,1:C,M); axis tight xy; colorbar
    title(sprintf('Band %d',k));
end
xlabel('Time (s)'); ylabel('Channel');
end
