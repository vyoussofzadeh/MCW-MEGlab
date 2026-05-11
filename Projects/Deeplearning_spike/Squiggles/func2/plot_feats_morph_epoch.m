function plot_feats_morph_epoch(X, fs, chIdx)
% X: C×T, already prepped; chIdx optional
M = feats_morph_time(X,fs);
t=(0:size(X,2)-1)/fs; C=size(X,1);
figure('Color','w'); tiledlayout(5,1,'TileSpacing','compact','Padding','compact')
nexttile; imagesc(t,1:C,X);     axis tight xy; colorbar; title('raw'); ylabel('chan')
nexttile; imagesc(t,1:C,M.d1);  axis tight xy; colorbar; title('d1')
nexttile; imagesc(t,1:C,M.d2);  axis tight xy; colorbar; title('d2')
nexttile; imagesc(t,1:C,M.env); axis tight xy; colorbar; title('env'); xlabel('time (s)')
nexttile; bar(M.zcr); xlim([0 C+1]); grid on; ylabel('ZCR'); title('zero-crossings')

if nargin>=3 && ~isempty(chIdx)
    figure('Color','w'); 
    plot(t,X(chIdx,:),'k'); hold on; plot(t,M.d1(chIdx,:),'b'); plot(t,M.d2(chIdx,:),'r');
    plot(t, rescale(M.env(chIdx,:), min(X(chIdx,:)), max(X(chIdx,:))), 'g');
    legend('raw','d1','d2','env(rescaled)'); grid on
    title(sprintf('Channel %d overlay', chIdx));
end
end
