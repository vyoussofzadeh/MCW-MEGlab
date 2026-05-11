function plot_mean_feature_maps(Xcells, Y, fs, names)
D = size(Xcells{1},3); C = size(Xcells{1},1); T = size(Xcells{1},2);
mS = zeros(C,T,D); nS = 0; mN = zeros(C,T,D); nN = 0;
for i=1:numel(Xcells)
    if Y(i)=="Spike",   mS = mS + Xcells{i}; nS=nS+1;
    else,               mN = mN + Xcells{i}; nN=nN+1;
    end
end
mS = mS/max(nS,1); mN = mN/max(nN,1);

t = (0:T-1)/fs;
figure('Color','w'); tiledlayout(2,D,'TileSpacing','compact','Padding','compact');
for k=1:D
    nexttile; imagesc(t,1:C,mS(:,:,k)); axis tight xy; title(['Spike-' names{k}]); colormap(parula); colorbar
    nexttile; imagesc(t,1:C,mN(:,:,k)); axis tight xy; title(['NoSpike-' names{k}]); colorbar
end
end

% Example:
% plot_mean_feature_maps(Xtr2D, YTrain, fs, {'raw','d1','d2','env'});
