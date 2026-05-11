% function [X2D, names] = make_stack_2D_spectral_moving(X, fs, bands, winSec, hopSec)
% % MAKE_STACK_2D_SPECTRAL_MOVING  Build C×T×D stack: raw + moving bandpower maps
% % X   : C×T epoch (already filtered/scaled/norm/padded)
% % fs  : sampling rate (Hz)
% % bands : K×2 [low high] Hz (default d..?)
% % winSec, hopSec : moving window/hop in seconds (defaults 0.25 / 0.05)
% %
% % Returns:
% %   X2D  : C×T×D, where D = 1 + K  (raw + one map per band)
% %   names: 1×D cell array of strings, e.g., {'raw','delta','theta',...}
% 
%     if nargin<3 || isempty(bands)
%         bands = [0.5 4; 4 8; 8 13; 13 30; 30 70];  % d ? a ß ?
%     end
%     if nargin<4 || isempty(winSec), winSec = 0.25; end
%     if nargin<5 || isempty(hopSec), hopSec = 0.05; end
% 
%     % Compute moving bandpower maps (C×T for each band)
%     B = feats_spectral_moving(X, fs, bands, winSec, hopSec);
% 
%     % Stack dynamically: raw first, then each band map in order
%     K = size(bands,1);
%     maps = cell(1, K);
%     for k = 1:K
%         fn = sprintf('bp%d', k);
%         if ~isfield(B, fn)
%             error('make_stack_2D_spectral_moving: missing field %s in B.', fn);
%         end
%         maps{k} = B.(fn);
%     end
%     X2D = single(cat(3, X, maps{:}));
% 
%     % Build names safely
%     defaultNames = {'delta','theta','alpha','beta','gamma','high'};
%     bandNames = defaultNames(1:min(K, numel(defaultNames)));
%     % If K > numel(defaultNames), auto-name the extras:
%     if K > numel(defaultNames)
%         for kk = numel(defaultNames)+1 : K
%             bandNames{kk} = sprintf('band%d', kk); %#ok<AGROW>
%         end
%     end
%     names = [{'raw'}, bandNames];
% end


function [X2D, names] = make_stack_2D_spectral_moving(X, fs, bands, winSec, hopSec)
if nargin<3 || isempty(bands)
    bands = [0.5 4; 4 8; 8 13; 13 30; 30 70];
end
if nargin<4, winSec = 0.25; end
if nargin<5, hopSec = 0.05; end

B = feats_spectral_moving(X, fs, bands, winSec, hopSec);
K = size(bands,1);

maps = cell(1,K);
for k=1:K, maps{k} = B.(sprintf('bp%d',k)); end
X2D = single(cat(3, X, maps{:}));          % C×T×(1+K)

defaultNames = {'delta','theta','alpha','beta','gamma','high'};
bandNames = defaultNames(1:min(K, numel(defaultNames)));
if K > numel(defaultNames)
    for kk=numel(defaultNames)+1:K, bandNames{kk} = sprintf('band%d',kk); end
end
names = [{'raw'}, bandNames];
end
