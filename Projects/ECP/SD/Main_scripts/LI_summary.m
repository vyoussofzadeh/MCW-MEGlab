% -----------------------------------------------------------
% Concordance bar chart: fixed- vs subject-specific windows
% -----------------------------------------------------------

% Data (rows = LI methods; columns = four bar categories)
% Columns: 1 = Animal-Baseline | Fixed
%          2 = Animal-Baseline | Subject-specific
%          3 = Animal-Symbol   | Fixed
%          4 = Animal-Symbol   | Subject-specific
% M = [ 74.1  78.1  88.8  90.2 ;   % Magnitude
%       73.2  77.6  87.5  89.3 ;   % Counting
%       71.8  76.8  86.1  88.6 ];  % Bootstrap
 
  
m1 = [67.6   71.8
   71.8   70.4
   71.8   71.8];

m2 = [85.1   89.1
   87.8   89.1
   89.1   90.5];

M = [m1, m2];
   

methods   = {'Magnitude','Counting','Bootstrap'};

barLabels = { 'A-vs-B | Fixed' ...
    'A-vs-B | Peak'  ...
    'A-vs-S | Fixed' ...
    'A-vs-S | Peak'  };

% Colour-blind-safe palette (4 columns)
cmap = [ 0.96 0.60 0.15 ;   % orange
         0.85 0.40 0.18 ;   % dark orange
         0.40 0.76 0.65 ;   % teal
         0.23 0.50 0.70 ];  % blue-teal

figure('Units','pixels','Position',[100 100 650 320])
b = bar(M,'grouped'); hold on
for k = 1:numel(b)
    b(k).FaceColor = cmap(k,:);
    b(k).EdgeColor = 'none';
end

% Y-axis
ylabel('MEG - fMRI LI concordance (%)','FontSize',11)
ylim([60 100]); yticks(60:10:100); grid on
ax = gca; ax.GridLineStyle=':'; ax.YGrid='on';
ax.XTick = 1:numel(methods); ax.XTickLabel = methods;
ax.FontSize = 10; ax.Box = 'off';

% Legend
legend(barLabels,'Location','northoutside','Orientation','horizontal',...
       'FontSize',9,'Box','off')

% Value labels inside bars
for i = 1:size(M,1)
    for j = 1:size(M,2)
        x = b(j).XData(i) + b(j).XOffset;
        y = M(i,j);
        text(x, y/2, sprintf('%.1f%%',y), ...
             'HorizontalAlignment','center', ...
             'VerticalAlignment','middle', ...
             'Color','w','FontSize',8,'FontWeight','bold')
    end
end
set(gca,'color','none');


% title({'Fixed intervals vs. subject-specific peak windows'; ...
%        'Effect of control-task contrast on concordance'}, ...
%        'FontSize',12,'FontWeight','normal')
