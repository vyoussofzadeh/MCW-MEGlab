% plotSNR(final_combined_snr.mean_tSSSvsRaw, string(1:length(final_combined_snr.SubjectID)), ...
%     'Median SNR Across Trials, tSSS-vs-Raw', [0.4922 0.0039 0.9063], [0.5, 0.5, 0.5], 100);
% doPlotExport(plot_option, save_dir, 'MEGnetvstSSS', 'svg');
% 
% plotSNR(final_combined_snr.mean_MEGnetvstSSS, string(1:length(final_combined_snr.SubjectID)), ...
%     'Median SNR Across Trials, tSSS-vs-Raw', [0.9922 0.8672 0.0039], [0.5, 0.5, 0.5], 100);
% doPlotExport(plot_option, save_dir, 'tSSSvsRaw', 'svg');

DataArray = [final_combined_snr.nSNR_tSSSvsRaw, final_combined_snr.nSNR_MEGnetvstSSS];
Colors = [0.4922 0.0039 0.9063; 0.9922 0.8672 0.0039];

figure; hold on;
[xPositions, yPositions, ~, ~] = UnivarScatter(DataArray,'Label',{'Raw vs tSSS','tSSS vs MEGnet'},'MarkerFaceColor',Colors, 'StdColor',[0.5, 0.5, 0.5],'SEMColor', [0.7,0.7,0.7],'PointSize',103);
ylabel('SNR (dB) Median','FontSize', 16);
xlabel('Subjects','FontSize', 16);
set(gca,'color','none');
title('SNR (dB) with MEG','fontsize',16)
set(gca,'FontName','HelveticaNeueLT Std Lt');
ylim([-45 20]);

f = [xPositions, yPositions];
for j=1:length(f)
    line([f(j,1),f(j,2)],[f(j,3),f(j,4)], 'LineWidth', 1); % Adjust line width if necessary
    text(xPositions(j,1)+0.05,yPositions(j,1), num2str(j))
    text(xPositions(j,2)+0.05,yPositions(j,2), num2str(j))
end
xlim([0.5 2.5]);
doPlotExport(plot_option, save_dir, 'Scatter_noiseSNR', 'svg');

%
% Number of subjects (assuming both vectors have the same length)
nSubjects = length(final_combined_snr.SubjectID);

% Open a new figure
figure;
hold on;  % allows multiple scatter plots on the same axes

% First scatter: tSSS vs Raw
scatter(1:nSubjects, final_combined_snr.nSNR_tSSSvsRaw, ...
    'SizeData', 100, ...                      % Marker size (in points^2)
    'MarkerEdgeColor', [0.5, 0.5, 0.5], ...    % Border color (gray)
    'MarkerFaceColor', [0.4922, 0.0039, 0.9063]);  % Fill color (purple-ish)

% Second scatter: MEGnet vs tSSS
scatter(1:nSubjects, final_combined_snr.nSNR_MEGnetvstSSS, ...
    'SizeData', 100, ...
    'MarkerEdgeColor', [0.5, 0.5, 0.5], ...
    'MarkerFaceColor', [0.9922, 0.8672, 0.0039]); % Fill color (yellow-ish)

% Set up axis labels, ticks, etc.
xticks(1:nSubjects);
xticklabels(string(1:nSubjects));
xlim([0.5, nSubjects + 0.5]);
set(gca,'FontName','HelveticaNeueLT Std Lt');

% Add a legend so you can tell which points are which
legend('Raw vs tSSS','tSSS vs MEGnet','Location','best');

% Some final cleanup
set(gca, 'FontSize', 8);
set(gca, 'color', 'none'); 
set(gcf, 'Position', [700, 500, 1500, 300]);
title('Noise SNR: tSSS-vs-Raw vs MEGnet-vs-SSS', 'FontSize', 16);
xlabel('Subjects','FontSize', 16);
ylabel('SNR (dB)','FontSize', 16);

% Optional: export the figure
doPlotExport(plot_option, save_dir, 'noiseSNR', 'svg');
