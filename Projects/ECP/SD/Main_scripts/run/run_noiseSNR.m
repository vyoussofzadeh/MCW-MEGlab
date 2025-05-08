%% MEG noise-SNR
% Load Data and Compute med_snr_mag / med_snr_grad
cd('/data/MEG/Research/MEGIN project/Scripts/SNR/SNR_Data')
[all_snr_grad, all_snr_mag] = ecpfunc_read_noiseSNR();

snr_MEGnetvstSSS  = [all_snr_mag(1).snr_value,  all_snr_grad(1).snr_value];
snr_tSSSvsRaw = [all_snr_mag(2).snr_value, all_snr_grad(2).snr_value];

%%
% Reshape the 3D array into a 2D array with dimensions 10x30
snr_MEGnetvstSSS1 = reshape(snr_MEGnetvstSSS, size(snr_MEGnetvstSSS, 1), []);
snr_MEGnetvstSSS2 = median(snr_MEGnetvstSSS1,2);

snr_tSSSvsRaw1 = reshape(snr_tSSSvsRaw, size(snr_tSSSvsRaw, 1), []);
snr_tSSSvsRaw2 = median(snr_tSSSvsRaw1,2);

%%
% med_snr_tSSSvsRaw  = squeeze(median(snr_MEGnetvstSSS, 2));
% med_snr_MEGnetvstSSS = squeeze(median(snr_tSSSvsRaw, 2));

% mean_tSSSvsRaw  = mean(med_snr_tSSSvsRaw,  2);
% mean_MEGnetvstSSS = mean(med_snr_MEGnetvstSSS, 2);

% Identify or retrieve subject IDs:
sub_snr_list = unique(all_snr_mag(1).sub);  % cell array of length N_subj
% if length(sub_snr_list) ~= nSubjects
%     error('Mismatch: sub_snr_list has %d entries, but mean_snr has %d.', ...
%         length(sub_snr_list), nSubjects);
% end

plot_option = 0;
if plot_option ==1
    
    plotSNR(snr_tSSSvsRaw2, sub_snr_list, 'Mean SNR Across Subjects, tSSSvsRaw');
    doPlotExport(plot_option, save_dir, 'tSSSvsRaw', 'svg');

    plotSNR(snr_MEGnetvstSSS2, sub_snr_list, 'Mean SNR Across Subjects, MEGnetvstSSS');
    doPlotExport(plot_option, save_dir, 'MEGnetvstSSS', 'svg');
    
    %% all data put in one matrix_1 column_tsss and 2_column_megnet
    DataArray = [snr_tSSSvsRaw2, snr_MEGnetvstSSS2];
    
    %% Customized Scatter Plot with Connected Data Points Across Groups
    % addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/other')
    Colors = [0.4922 0.0039 0.9063; 0.9922 0.8672 0.0039];
    
    %%
    figure;
    hold on;
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

    
end

%%
% Build a table with 4 columns:
T_snr = table(sub_snr_list(:), snr_tSSSvsRaw2, snr_MEGnetvstSSS2, ...
    'VariableNames', {'SubjectID','nSNR_tSSSvsRaw','nSNR_MEGnetvstSSS'});

% 1) Convert T_snr.SubjectID to uppercase
if iscell(T_snr.SubjectID)
    % If it's a cell array of char
    T_snr.SubjectID = upper(string(T_snr.SubjectID));
else
    % If it's already a string array
    T_snr.SubjectID = upper(T_snr.SubjectID);
end

