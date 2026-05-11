function plot_stacked_vertical(X, fs, chs)
% X: [C x T], fs in Hz, chs (optional): indices to plot (default all)
if nargin<3 || isempty(chs), chs = 1:size(X,1); end
X = X(chs,:); [C,T] = size(X);
t = (0:T-1)/fs;

% normalize each row and add vertical offsets
x = X ./ (std(X,0,2)+eps);                 % z-scored per channel
off = (0:C-1)';                             % 0,1,2,...
Y = x + off;                                % shift each row up

figure('Color','w'); plot(t, Y', 'k'); grid on; xlim([t(1) t(end)])
yticks(off); yticklabels(string(chs));
ylabel('Channel'); xlabel('Time (s)');
title(sprintf('Stacked traces (%d channels)', numel(chs)));
end
