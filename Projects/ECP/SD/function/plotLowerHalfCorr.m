function plotLowerHalfCorr(T, numericVars, numericVarsShort, tit)
% Extract data
X = T{:, numericVars};

% Suppose you have the correlation matrix R and p-values p
[R, p] = corr(X, 'Rows','pairwise');

% Create a lower-triangle mask (including diagonal)
nVars = size(R,1);
lowerMask = tril(true(nVars));   % logical mask for i>=j

% Make a copy of R (and set the upper triangle to NaN)
Rlower = R;
Rlower(~lowerMask) = NaN;

% Create a figure
figure('Color','w','Position',[200,200,900,700]);

% Plot using imagesc
hImg = imagesc(Rlower, [-1,1]);
axis square
colormap jet
colorbar

% --- Make the NaN pixels transparent ---
set(hImg, 'AlphaData', ~isnan(Rlower));
set(gca, 'Color','w');  % the axes background is white

title(tit);

% Adjust ticks/labels
set(gca,...
    'XTick',1:nVars, 'XTickLabel',numericVarsShort,...
    'YTick',1:nVars, 'YTickLabel',numericVarsShort,...
    'XTickLabelRotation',45,...
    'TickLabelInterpreter','none',...
    'FontSize',9);

% Label only where you want (e.g., significant cells, or all in lower half)
hold on
alphaLevel = 0.05;
for i = 1:nVars
    for j = 1:i
        if p(i,j) < alphaLevel
            rVal = R(i,j);
            labelStr = sprintf('%.2f*', rVal);
            % Decide text color based on correlation
            txtColor = 'k';
            if abs(rVal) > 0.5
                txtColor = 'w';
            end
            text(j, i, labelStr, 'Color',txtColor, ...
                'HorizontalAlignment','center','FontSize',8);
        end
    end
end

%     % Draw group boundaries
%     for b = 1:numel(groupBoundaries)
%         loc = groupBoundaries(b) + 0.5;
%         xline(loc, 'Color','k', 'LineWidth',1.5);
%         yline(loc, 'Color','k', 'LineWidth',1.5);
%     end

hold off

end

% --- Simple custom diverging colormap function ---
function cmap = bwr_colormap(n)
if nargin<1, n = 256; end
half = floor(n/2);
% Blue to white
bluetowhite = [linspace(0,1,half)', linspace(0,1,half)', ones(1,half)'];
% White to red
whitetored  = [ones(1, n-half)', linspace(1,0,n-half)', linspace(1,0,n-half)'];
cmap = [bluetowhite; whitetored];
end
