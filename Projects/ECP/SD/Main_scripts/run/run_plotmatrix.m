% ========================= Example Plotmatrix Script =========================
%
% This script assumes you have:
%   1) A table named T in the workspace
%   2) numericVars:    a cell array of variable names (strings) in T
%   3) numericVarsLabel: a parallel cell array of labels to appear on the axes

%% 1) User-defined variable inputs
numericVars = { ...
    'optMEG_LI','fMRI_LI','rSNR','Animal_RT','Symbol_RT', ...
    'Animal_ACC','Symbol_ACC','AEDCount','EHQ','CP_freq','FSIQ'};
numericVarsLabel = { ...
    'MEG-LI','fMRI-LI','rSNR','Anim-RT','Symb-RT', ...
    'Anim-ACC','Symb-ACC','AEDCount','EHQ','CP-freq','FSIQ'};

% Extract numeric data from table T
X = T{:, numericVars};

%% 2) Generate the scatter matrix
% h: MxM array of graphics handles (off-diagonal => scatter, diagonal => hist)
% AX: MxM array of axes handles
[h, AX] = plotmatrix(X);

% Number of variables
nVars = size(X, 2);

%% 3) Label the axes
%   - The bottom row: x-labels
%   - The left column: y-labels
for i = 1:nVars
    AX(nVars, i).XLabel.String = numericVarsLabel{i};  % x-label
    AX(i, 1).YLabel.String = numericVarsLabel{i};      % y-label
end

%% 4) Make off-diagonal scatter markers unfilled (hollow circles)
%   - The diagonal subplots are histograms, so skip (i == j).
for row = 1:nVars
    for col = 1:nVars
        if row ~= col
            set(h(row,col), ...
                'Marker','o', ...            % circle marker
                'MarkerFaceColor','none', ...% hollow center
                'MarkerEdgeColor','k', ...   % black outline
                'LineStyle','none', ...      % ensures scatter, no connecting line
                'MarkerSize',2);            % optional marker size
        end
    end
end

%% 5) Draw x=0 and y=0 dashed lines in off-diagonal subplots
for row = 1:nVars
    for col = 1:nVars
        if row ~= col
            ax = AX(row,col);
            hold(ax, 'on');
            xline(ax, 0,'k--','LineWidth',0.8,'HandleVisibility','off');
            yline(ax, 0,'k--','LineWidth',0.8,'HandleVisibility','off');
        end
    end
end

%% 6) Reduce font sizes for tick labels and axis labels
for row = 1:nVars
    for col = 1:nVars
        axHandle = AX(row, col);
        axHandle.FontSize = 8;          % Tick label font size
        axHandle.XLabel.FontSize = 8;   % X-axis label font size
        axHandle.YLabel.FontSize = 8;   % Y-axis label font size
    end
end

%% 7) Final figure styling
set(gcf, 'Name','Pairwise Scatter Matrix','NumberTitle','off', 'Color','w');
