function M = bandmaps_from_TFR(Pow, f, bands)
% Pow: C×F×T; f: F×1; bands: K×2 Hz ? fields M.bp1..bpK (C×T)
if nargin<3, bands = [0.5 4; 4 8; 8 13; 13 30; 30 70]; end
C=size(Pow,1); T=size(Pow,3); K=size(bands,1);
for k=1:K
    m  = f>=bands(k,1) & f<bands(k,2);
    M.(sprintf('bp%d',k)) = squeeze(sum(Pow(:,m,:),2));   % integrate power across freqs ? C×T
end
end
