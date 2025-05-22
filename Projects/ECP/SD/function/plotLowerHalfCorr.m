function plotLowerHalfCorr(T, numericVars, numericVarsShort, tit)
% Extract data
X = T{:, numericVars};

% Suppose you have the correlation matrix R and p-values p
[R, p] = corr(X, 'Rows','pairwise');


% ------------------------------------------------------------------
% 2. BenjaminiHochberg FDR on the lower-triangle p-values
nVars   = size(R,1);
triMask = tril(true(nVars),-1);          % lower-triangle w/o diagonal
pVec    = p(triMask);                    % vector of raw p-values

qVec = mafdr(pVec,'BHFDR',true);         % q-values (needs Statistics TBX)
% --- fallback if mafdr is missing ---
% [ps,idx] = sort(pVec); m=numel(ps);
% qtmp = ps.*m./(1:m)'; qtmp = min(cummin(flipud(qtmp)),1);
% qVec(idx) = qtmp;

q           = nan(nVars);                % same size as R
q(triMask)  = qVec;

%%

% Create a lower-triangle mask (including diagonal)
nVars = size(R,1);
% lowerMask = tril(true(nVars));   % logical mask for i>=j
lowerMask = tril(true(nVars),-1);   % <- use 1 offset to exclude diag

% Make a copy of R (and set the upper triangle to NaN)
% Rlower = R;
% Rlower(~lowerMask) = NaN;

% Apply the mask
Rlower          = NaN(nVars);   % start with all NaNs
Rlower(lowerMask) = R(lowerMask);

% Create a figure
figure('Color','w','Position',[200,200,900,700]);

% Plot using imagesc
hImg = imagesc(Rlower, [-1,1]);
axis square
colormap jet
colormap(bwr_colormap(255));                 % symmetric diverging map
% colorbar

% --- Make the NaN pixels transparent ---
set(hImg, 'AlphaData', ~isnan(Rlower));
set(gca, 'Color','w');  % the axes background is white

% ----------------- move main heat-map up a bit --------------------------
ax      = gca;
axPos   = ax.Position;        % [left bottom width height]  (normalized)
gap     = 0.08;               % 8 % of figure height
axPos(2)= axPos(2) + gap;     % raise the axes
axPos(4)= axPos(4) - gap;     % keep square shape (shrink height)
ax.Position = axPos;

% ----------------- resize + drop the colour-bar -------------------------
cb      = colorbar('southoutside');
cb.Ticks = -1:0.5:1;
cb.Label.String = 'Pearson r';

cbPos   = cb.Position;        % [left bottom width height]
cbPos(2)= cbPos(2) - 0.18;    % push it 4 % lower (avoids overlap)
cbPos(3)= cbPos(3)*0.60;      % 60 % as wide  ? slimmer bar
cbPos(1)= axPos(1) + (axPos(3)-cbPos(3))/2;   % centre under plot
cb.Position = cbPos;

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
for i = 2:nVars
    for j = 1:i-1
        if q(i,j) < alphaLevel        % ? use FDR-corrected q % p(i,j) < alphaLevel
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

% Adjust ticks/labels
set(gca,...
    'XTick',1:nVars,'XTickLabel',numericVarsShort,...
    'YTick',1:nVars,'YTickLabel',numericVarsShort,...
    'XTickLabelRotation',45,...
    'TickLabelInterpreter','none',...
    'FontSize',9);

% --- hide the axes box & ticks but keep labels --------------------------
ax = gca;
ax.Box        = 'off';      % removes the box border
ax.TickLength = [0 0];      % hides tick marks
ax.XRuler.Axle.Visible = 'off';   % hide X-axis spine
ax.YRuler.Axle.Visible = 'off';   % hide Y-axis spine


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
