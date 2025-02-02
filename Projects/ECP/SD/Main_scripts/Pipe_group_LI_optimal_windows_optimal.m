
%% ECP Semantic decision task dataset, Medical College of Wisconsin

% Script: BS Process (Laterality analysis)
% Project: ECP_SD
% Writtern by: Vahab Youssof Zadeh
% Update: 05/21/2024

clear; clc, close('all'); warning off,

%% Paths
restoredefaultpath
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/run')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/function')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/tools/helpful_tools/daviolinplot/daboxplot')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/Main_scripts/run')
Run_setpath

addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/other/');


%%
MEG_thre = 10; % MEG threshold
fMRI_thre = 10; % fMRI threshold

figtype = 'svg';

%%
LI_analysis_label = {'DICS_indirect','DICS_directcontrast','LCMV_anim_vs_Symb','-','DICS_anim', 'DICS_contrast_prestim', 'dSPM_contrast','DICS_directcontrast_localthresh'};

for i = 1:length(LI_analysis_label)
    disp([num2str(i) ') ' LI_analysis_label{i}]);
end

LI_analysis = input('');

%% plotting option
disp('1: plot in browser (and save in .svg), 2: plot in matlab figure (not saved)')
plot_option = input('');

%%
LI_method_label = {'Magnitude', 'Counting','Bootstrapping'};

disp('1: Magnitude')
disp('2: Counting')
disp('3: Bootstrapping')
% LI_method = input(':');

%%
data_save_dir = '/data/MEG/Research/ECP/Semantic_Decision/Results_prestim';
cd(data_save_dir)

%%
clc
save_dir = fullfile(data_save_dir,LI_analysis_label{LI_analysis}, 'compare_LIs');
checkOrCreateDir(save_dir)
cd(save_dir)

switch LI_analysis
    case 3
        LI_pt_run = []; rSNR_run = [];
         % Here we have two runs to process
        for trun = 1:2
            % Create a label for each run
            runLabel = sprintf('run%d', trun);

            for i = 1:length(LI_method_label)
                currentMethod = LI_method_label{i};

                % --- Load the data ---
                loadedData = load(fullfile(...
                    data_save_dir, ...
                    LI_analysis_label{LI_analysis}, ...
                    [currentMethod, '_run', num2str(trun)], ...
                    'LI_Patn' ...
                ));

                % --- Store in a sub-struct for each run ---
                LI_pt_run.(runLabel).(currentMethod) = loadedData;

                % Decide if we store pow_sub vs. count_sub
                if i == 1
                    % If i=1 -> store 'pow_sub'
                    rSNR_run.(runLabel).(currentMethod) = loadedData.pow_sub;
                else
                    % If i>1 -> store 'count_sub'
%                     rSNR.(runLabel).(currentMethod) = loadedData.count_sub;
                end
            end
        end
        
                % --- CHOOSE the run for the rest of the analysis ---
        % For example, if chosenRun = 1, analyze run1 only
        disp('choose runs (1 and 2)');
        chosenRun = input('');
        selectedRunLabel = sprintf('run%d', chosenRun);

        % Now pass only the selected run's data to the rest of your analysis
        LI_pt = LI_pt_run.(selectedRunLabel);
        rSNR  = rSNR_run.(selectedRunLabel);

        % Next, do your analysis with LI_pt_selected / rSNR_selected.
        % Example:
        disp(['Using run ', num2str(chosenRun), ' for analysis...']);
        % e.g. do_something(LI_pt_selected, rSNR_selected);
    otherwise
        LI_pt = []; rSNR = [];
        for i=1:length(LI_method_label)
            LI_pt.(LI_method_label{i}) = load(fullfile(data_save_dir,LI_analysis_label{LI_analysis},LI_method_label{i},'LI_Patn'));
            if i == 1
                rSNR.(LI_method_label{i}) = LI_pt.(LI_method_label{i}).pow_sub;
            else
                rSNR.(LI_method_label{i}) = LI_pt.(LI_method_label{i}).count_sub;
            end
        end
end

%% Run Brainstorm
bs_path = '/opt/matlab_toolboxes/Brainstorm/Brainstorm3_2022/brainstorm3';

BS_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/';
protocol = fullfile(BS_dir, 'data_full/protocol.mat');

%%
Run_load_surface_template

%% Read fMRI lats
fmri_LIs = ecpfunc_read_fmri_lat();

%%
datadir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data';
BS_data_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/data_full';

cfg = [];
cfg.protocol = protocol;
cfg.datadir = datadir;
cfg.BS_data_dir = BS_data_dir;

switch LI_analysis
    case 1
        cfg.datatag = 'wDICS_baseline_18_4_50ms';
        S_data = ecpfunc_read_sourcemaps_dics(cfg);
    case {2,8}
        cfg.datatag = 'wDICS_18_4_50ms';
        S_data = ecpfunc_read_sourcemaps_dics_contrast(cfg);
    case 3
        %         cfg.datamask = fullfile('./Group_analysis/LCMV/results_average*.mat');
        %         S_data = ecpfunc_read_sourcemaps(cfg);
        protocol = fullfile(BS_dir, 'data_full/protocol.mat');
        cfg.protocol = protocol;
        BS_data_dir = fullfile(BS_dir,'data');
        cfg.BS_data_dir = BS_data_dir;
        
        cfg.datamask =  './Group_analysis/ec*/*abs_ssmooth.mat';
        S_data = ecpfunc_read_sourcemaps2(cfg);
        S_data_sep = ecpfunc_separate_runs(S_data);
    case 5
        protocol = fullfile(BS_dir, 'data/protocol.mat');
        BS_data_dir = fullfile(BS_dir,'data');
        cfg.protocol = protocol;
        cfg.BS_data_dir = BS_data_dir;
        cfg.datatag = 'wDICS_baseline_18_4';
        S_data = ecpfunc_read_sourcemaps_dics_anim(cfg);
    case 6
        cfg.datatag = 'wDICS_contrast_18_4_50ms';
        S_data = ecpfunc_read_sourcemaps_dics_contrast(cfg);
    case 7
        cfg.datatag = 'dSPM_contrast';
        S_data = ecpfunc_read_sourcemaps_contrast(cfg);
end

%% Subject demog details
switch LI_analysis
    case {1,5}
        cfg = []; cfg.subjs_3 = S_data.subjs_3; cfg.subjs_2 = S_data.subjs_2;
        cfg.sFiles_3 = S_data.sFiles_3; cfg.sFiles_2 = S_data.sFiles_2;
        sub_demog_data = ecpfunc_read_sub_demog(cfg);
    case 3
        cfg = [];
        cfg.subjs_3  = S_data_sep.subjs_3_run1;
        cfg.subjs_2  = S_data_sep.subjs_2_run1;
        cfg.sFiles_3 = S_data_sep.sFiles_3_run1;
        cfg.sFiles_2 = S_data_sep.sFiles_2_run1;
        sub_demog_data_run1 = ecpfunc_read_sub_demog(cfg);
        
        cfg = [];
        cfg.subjs_3  = S_data_sep.subjs_3_run2;
        cfg.subjs_2  = S_data_sep.subjs_2_run2;
        cfg.sFiles_3 = S_data_sep.sFiles_3_run2;
        cfg.sFiles_2 = S_data_sep.sFiles_2_run2;
        sub_demog_data_run2 = ecpfunc_read_sub_demog(cfg);
    case {2,4,6,8}
        %         clc
        cfg = []; cfg.subjs = S_data.subjs;
        cfg.sFiles = S_data.sFiles_32;
        sub_demog_data = ecpfunc_read_sub_demog_contrast(cfg);
end

%% TLE side (PT only)
patn_neuropsych_data = ecpfunc_read_patn_neuropsych();

%%
switch LI_analysis
    case {1,5}
        cfg = [];
        cfg.sub_demog_data = sub_demog_data;
        cfg.patn_neuropsych_data = patn_neuropsych_data;
        sub_TLE_sub_data = ecpfunc_read_sub_TLE_sub(cfg);
    case 3
        sub_demog_data_runs = {sub_demog_data_run1, sub_demog_data_run2};
        
        cfg = [];
        cfg.sub_demog_data = sub_demog_data_runs{chosenRun};
        cfg.patn_neuropsych_data = patn_neuropsych_data;
        sub_TLE_sub_data = ecpfunc_read_sub_TLE_sub(cfg);
        
    case {2,4,6,8}
        %         clc
        cfg = [];
        cfg.sub_demog_data = sub_demog_data;
        cfg.patn_neuropsych_data = patn_neuropsych_data;
        sub_TLE_sub_data = ecpfunc_read_sub_TLE_sub_contrast(cfg);
end

%% Inter-subject (group) averaging,
switch LI_analysis
    case {1,5}
        disp('1: Anim, Ctrl')
        disp('2: Anim, Patn')
        disp('3: Symbol, Ctrl')
        disp('4: Symbol, Patn')
        cfg = [];
        cfg.sub_demog_data = sub_demog_data;
        cfg.select_data = 1; S_data_anim_hc = ecpfunc_select_data(cfg);
        cfg.select_data = 2; S_data_anim_pt = ecpfunc_select_data(cfg);
        cfg.select_data = 3; S_data_symb_hc = ecpfunc_select_data(cfg);
        cfg.select_data = 4; S_data_symb_pt = ecpfunc_select_data(cfg);
        [LI_pt_ID,~,~] = intersect(S_data_anim_pt.sFiles_subid, S_data_symb_pt.sFiles_subid);
        [LI_hc_ID,~,~] = intersect(S_data_anim_hc.sFiles_subid, S_data_symb_hc.sFiles_subid);
    case 3
        disp('1: Anim, Ctrl')
        disp('2: Anim, Patn')
        disp('3: Symbol, Ctrl')
        disp('4: Symbol, Patn')
        cfg = [];
        cfg.sub_demog_data = sub_demog_data_runs{chosenRun};
        cfg.select_data = 1; S_data_anim_hc = ecpfunc_select_data(cfg);
        cfg.select_data = 2; S_data_anim_pt = ecpfunc_select_data(cfg);
        cfg.select_data = 3; S_data_symb_hc = ecpfunc_select_data(cfg);
        cfg.select_data = 4; S_data_symb_pt = ecpfunc_select_data(cfg);
        [LI_pt_ID,~,~] = intersect(S_data_anim_pt.sFiles_subid, S_data_symb_pt.sFiles_subid);
        [LI_hc_ID,~,~] = intersect(S_data_anim_hc.sFiles_subid, S_data_symb_hc.sFiles_subid);
        
    case {2,4,8}
        disp('1: Ctrl')
        disp('2: Patn')
        cfg = [];
        cfg.sub_demog_data = sub_demog_data;
        cfg.select_data = 1; S_data_hc = ecpfunc_select_contrast_data(cfg);
        cfg.select_data = 2; S_data_pt = ecpfunc_select_contrast_data(cfg);
        LI_pt_ID = S_data_pt.sFiles_subid;
end

%%
cfg = []; cfg.strt = -0.5; cfg.spt = 2; cfg.overlap = 0.05; cfg.linterval = 0.3;
wi  = do_time_intervals(cfg);

%% HCP atlas
glass_atlas = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/Best/scout_mmp_in_mni_symmetrical_final_updated.mat';
src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/cortex_pial_low.fs';
glass_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/Glasser';

cfg = [];
cfg.src_fname = src_fname;
cfg.glass_dir = glass_dir;
cfg.glass_atlas = glass_atlas;
cfg.plotflag = 0;
Data_hcp_atlas = ecpfunc_hcp_atlas2(cfg);
net_sel_mutiple_label = Data_hcp_atlas.groups_labels';

net_sel_mutiple_label{1} = 'Ang';
net_sel_mutiple_label{2} = 'Front';
net_sel_mutiple_label{6} = 'Temp';
net_sel_mutiple_label{11} = 'Lat';

%%
network_sel = [1:3,6:11];
colr = distinguishable_colors(length(network_sel));

%%
patn_neuropsych_tle = ecpfunc_read_patn_neuropsych_tle();
TLESide = patn_neuropsych_tle.TLESide; SUBNO = patn_neuropsych_tle.SUBNO;

SUBNO_pt = [];
for i=1:length(LI_pt_ID)
    SUBNO_pt(i) = str2double(LI_pt_ID{i}(3:end));
end

[~,~,IB] = intersect(SUBNO_pt, SUBNO);
TLESide_sel = TLESide(IB);

TLE_left = find(TLESide_sel == 'Left');

%% MEG vs. fMRI lat analysis (PT)
[sub_MF_pt,IA,IB_megfmri] = intersect(LI_pt_ID, fmri_LIs.ID.language_Lateral');

LI_pt_val = LI_pt.(LI_method_label{1}).LI_sub;

fmri_LIs_val = (fmri_LIs.val.language_Lateral(IB_megfmri)); % Lateral regions was used.

% missing from fMRI
difference = setdiff(LI_pt_ID, sub_MF_pt');
disp('missing from fMRI')
disp(difference');

%% MEG LI vs fMRI LI (language_Lateral)
LI_pt_val_new = [];
rSNR_new = [];

for i=1:length(LI_method_label)
    LI_pt_val_new.(LI_method_label{i}) = LI_pt.(LI_method_label{i}).LI_sub(:, IA,:);
%     if i==1
%         rSNR_new.(LI_method_label{i}) = rSNR.(LI_method_label{i})(:, IA);
%     else
        rSNR_new.(LI_method_label{i}) = rSNR.(LI_method_label{i})(:, IA);
%     end
end

cfg = [];
cfg.thre = fMRI_thre; cfg.LI = fmri_LIs_val;
fmri_LIs_trn = do_ternary_classification2(cfg);
size(fmri_LIs_trn);

switch LI_analysis
    case 3
        Pt_ID  = LI_pt.Magnitude.setup.S_data{1}.sFiles_subid(IA)';
        Pt_sfiles  = LI_pt.Magnitude.setup.S_data{1}.sFiles_in(IA)';
    otherwise
        Pt_ID  = LI_pt.Magnitude.setup.S_data.sFiles_subid(IA)';
        Pt_sfiles  = LI_pt.Magnitude.setup.S_data.sFiles_in(IA)';
end

%%
% Loop through each element in the struct arrays
for i = 1:size(rSNR_new.Magnitude, 1)
    for j = 1:size(rSNR_new.Magnitude, 2)
        % Adjust Counting to make it a column vector
        if isrow(rSNR_new.Counting(i, j).left)
            rSNR_new.Counting(i, j).left = rSNR_new.Counting(i, j).left.';
        end
        
        % Adjust Counting for right if necessary
        if isrow(rSNR_new.Counting(i, j).right)
            rSNR_new.Counting(i, j).right = rSNR_new.Counting(i, j).right.';
        end
        
        % Calculate mean for Bootstrapping if it's not a single column vector
        if size(rSNR_new.Bootstrapping(i, j).left, 2) > 1
            rSNR_new.Bootstrapping(i, j).left = nanmean(rSNR_new.Bootstrapping(i, j).left, 2);
        end
        
        % Calculate mean for Bootstrapping right if necessary
        if size(rSNR_new.Bootstrapping(i, j).right, 2) > 1
            rSNR_new.Bootstrapping(i, j).right = nanmean(rSNR_new.Bootstrapping(i, j).right, 2);
        end
        
        % Ensure Magnitude is a column vector for both left and right
        if isrow(rSNR_new.Magnitude(i, j).left)
            rSNR_new.Magnitude(i, j).left = rSNR_new.Magnitude(i, j).left.';
        end
        if isrow(rSNR_new.Magnitude(i, j).right)
            rSNR_new.Magnitude(i, j).right = rSNR_new.Magnitude(i, j).right.';
        end
    end
end

%% Plot MEG LI for selected network ROIs
run_plot_MEGLIs

run_plot_MEGLIs_ver2

%%
run_plot_MEGsnr

%% Fixed interval analysis
%- Summerize LIs_fixedinterva
run_sumLIs_fixedinterval

%- plot_fixedinterval_corcon_LIsmethods
run_plot_fixedinterval_corcon_LIsmethods

%- plot_fixedinterval_corcon_rois1
run_plot_fixedinterval_corcon_rois1

%- plot_fixedinterval_corcon_rois2
% run_plot_fixedinterval_corcon_rois2

%- plot_fixedinterval_rois
run_plot_fixedinterval_rois

%- table_fixedinterval
run_table_fixedinterval

%- plot_fixedinterval
run_plot_fixedinterval

%% Dynamic interval analysis
plot_indiv_LI = 0; plot_rSNR = 0; plot_rSNR_LI = 0;

disp('1) res-snr and optimal bound selection')
disp('2) res-snr with fixed lower and upper bounds')
disp('3) res-snr with optimzation against fMRI LIs - only for inspection!')
disp('4) no optimzation against fMRI LIs (max fixed intervals)')
opt_sel = input('optimal method: ');

switch opt_sel
    case 1
        minlowerband = 0.35; maxUpperband = 1.5;
        %         minlowerband = 0.25; maxUpperband = 1.6;
        opt_method = 'rsnr_optbound'; %'rsnr_optbound_mean';
        %         opt_method = 'AUC'; %'rsnr_optbound_mean';
        run_optimalLIs_snr_dics_rois
        disp(['mean_conc:', num2str(mean(summaryTableDynamic_save.Concordance))])
    case 2
        lowerBound = 0.35; upperBound = 1.5;
        opt_method = 'rsnr'; run_optimalLIs_snr_dics
    case 3
        opt_method = 'rsnr'; run_optimalLIs_usingfMRI
    case 4
        opt_method = 'fixed';
        uniqueROIs = {'Ang', 'Front', 'Lat', 'Temp' };
        run_table_bestLIs_fixed_LI
        run_bestResultsTable_fixed
end

%%
% bestResultsTable_save.Best_LI_Method

% ecpfunc_blandAltmanPlot(bestResultsTable_save.optMEG_LI{1,1}, bestResultsTable.fMRI_LI{1,1});
% ecpfunc_blandAltmanPlot(bestResultsTable_save.optMEG_LI{2,1}, bestResultsTable.fMRI_LI{2,1});


%%
%- plot_compare_fixed_opt_methods
run_plot_compare_fixed_opt_methods

%- plot_compare_fixed_opt_rois
run_plot_compare_fixed_opt_rois

%- Plot_compare_fixed_opt
run_plot_compare_fixed_opt

%- plot_compare_fixed_opt_rois_summ
run_plot_compare_fixed_opt_rois_summ
run_plot_compare_fixed_opt_gsum

%% Pick the Best Results Out of 3 LI Methods for Discordant Analyses
run_table_bestLIs

% optimalMEG_LI = bestResultsTable.optMEG_LI{3};
% fMRI_LI = bestResultsTable.fMRI_LI{3};
% 
% % With labels:
% N = length(optimalMEG_LI);
% myLabels = arrayfun(@(x) sprintf('S%d', x), 1:N, 'UniformOutput', false);
% ecpfunc_blandAltmanPlot(optimalMEG_LI, fMRI_LI, myLabels);
% [concordanceBA, outlierIdx] = ecpfunc_blandAltmanConcordance(optimalMEG_LI, fMRI_LI);
% fprintf('Concordance = %.1f%%\n', 100*concordanceBA);
% fprintf('Outlier subjects: %s\n', num2str(outlierIdx'));

%% Optional 

close all
% Suppose your table is named "bestResultsTable" with columns:
%  - ROI
%  - Best_LI_Method
%  - MEG_LI  (each cell contains a [72 x 1] or [72 x T] numeric array)
%  - fMRI_LI (each cell contains a [72 x 1] or [72 x T] numeric array)
% etc.

for iRow = 1:height(bestResultsTable)
    
    % 1) Extract the ROI name and method (for labeling)
    currentROI    = bestResultsTable.ROI{iRow};
    currentMethod = bestResultsTable.Best_LI_Method{iRow};
    
    % 2) Extract the MEG_LI and fMRI_LI arrays
    %    Note: If these are [72 x 1], perfect; if they're bigger ([72 x 44]),
    %    you'll need to pick which column or reduce them somehow.
    megLI  = bestResultsTable.optMEG_LI{iRow};   % e.g., [72 x 1]
    fmriLI = bestResultsTable.fMRI_LI{iRow};  % e.g., [72 x 1]
    
    % 3) Generate subject labels (S1, S2, ...).
    %    Adjust if you have actual IDs in your data.
    nSubj  = size(megLI,1);  % e.g. 72
    subjectLabels = arrayfun(@(x) sprintf('S%d', x), 1:nSubj, 'UniformOutput', false);
    
    % 4) Create a new figure (optional)
    %     figure('Name', sprintf('BlandAltman: ROI=%s, Method=%s', currentROI, currentMethod));
    
    % 5) Call your BlandAltman plotting function.
    %    For example, if you want to label only outliers:
    ecpfunc_blandAltmanPlot(megLI, fmriLI, subjectLabels);
    [concordanceBA, outlierIdx] = ecpfunc_blandAltmanConcordance(megLI, fmriLI);
    disp([currentROI, '-', currentMethod])
    fprintf('Concordance = %.1f%%\n', 100*concordanceBA);
    fprintf('Outlier subjects: %s\n', num2str(outlierIdx'));
    
    % 6) (Optional) Adjust the plot title or save it
    title(sprintf('Bland-Altman Plot: %s (Method: %s)', currentROI, currentMethod));
    
    %    If you want to save automatically:
    %    saveas(gcf, sprintf('BlandAltman_%s_%s.png', currentROI, currentMethod));
end


%% Task and epil. measures
%- Response (reaction) time
run_responsereaction

% Task Performance
run_taskperformance

%- Epilepsy Metrics
T1_epil_measures = ecpfunc_read_epil_measures();   % Load additional measures

%% Combine all variables
combined1 = outerjoin(T_patn_MEGfMRI, accuracyResults_updt, 'Keys', 'SubjectID', 'MergeKeys', true);

% Merge the result with T1_epil_measures_tbl
final_combined = outerjoin(combined1, T1_epil_measures, 'LeftKeys', 'SubjectID', 'RightKeys', 'SubjectID', 'MergeKeys', true);
[commonSubs, idxFinal, idxPt] = intersect(final_combined.SubjectID, Pt_ID);
final_combined_updt       = final_combined(idxFinal,:);

%% Investigate Discordant Samples of Best Results and Obtain Corresponding MEG_LI and fMRI_LI
clc
close all

cfg = [];
cfg.roi_sel = 3; %lateral, n=3
cfg.wi = wi;
cfg.bounds = bounds;
cfg.plot_option  = 1;
cfg.save_dir = save_dir;
cfg.final_combined_updt = final_combined_updt;
cfg.T1_epil_measures = T1_epil_measures;
switch opt_sel
    case 1
        cfg.bestResultsTable = bestResultsTable;
    case 4
        cfg.bestResultsTable = bestResultsTable_fixed;
end
% ecpfunc_assess_discondances(cfg)
% cfg.discordSubs = [4, 5, 10, 21, 28, 36, 56];
% ecpfunc_assess_discondances(cfg)
ecpfunc_assess_gross_discondances(cfg)

%% Fisher analysis
run_fisheranalysis

%% Optional
% clc
ecpfunc_stackedBar_TLE_EHQ_IQ(bestResultsTable, T1_epil_measures)

%%
clc, close all

cfg = [];
cfg.bestResultsTable = bestResultsTable;
cfg.myCategorical = T1_epil_measures_upted.SymbolACCcat;
cfg.discordColumn = 'Gross_Discord_Subs';
cfg.categoryList  = {'Low','Mid','High'};
cfg.title = 'Symbol ACC';
cfg.doFisher = true;
cfg.binA = {'Low','Mid'};
cfg.binB = {'High'};
cfg.nSubjects = length(T1_epil_measures_upted.SubjectID);
[counts, pVals, ORvals, hFig] = ecpfunc_plotDiscordantStackedBar(cfg);
cfg.myCategorical = T1_epil_measures_upted.AnimalACCcat;
cfg.title = 'Animal ACC';
[counts, pVals, ORvals, hFig] = ecpfunc_plotDiscordantStackedBar(cfg);


cfg = [];
cfg.bestResultsTable = bestResultsTable;
cfg.myCategorical = T1_epil_measures_upted.AnimalRTcat;
cfg.discordColumn = 'Gross_Discord_Subs';
cfg.categoryList  = {'Fast','Moderate','Slow'};
cfg.doFisher = true;
cfg.binA = {'Moderate','Slow'};
cfg.binB = {'Fast'};
cfg.title = 'Animal RT';
cfg.nSubjects = length(T1_epil_measures_upted.SubjectID);
[counts, pVals, ORvals, hFig] = ecpfunc_plotDiscordantStackedBar(cfg);
cfg.myCategorical = T1_epil_measures_upted.SymbolRTcat;
cfg.title = 'Symbol RT';
[counts, pVals, ORvals, hFig] = ecpfunc_plotDiscordantStackedBar(cfg);


cfg = [];
cfg.bestResultsTable = bestResultsTable;
cfg.myCategorical = T1_epil_measures_upted.EHQcat;
cfg.discordColumn = 'Gross_Discord_Subs';
cfg.categoryList  = {'Left','Right','Ambi'};
cfg.doFisher = true;
cfg.binA = {'Ambi','Left'};
cfg.binB = {'Right'};
cfg.title = 'EHQ';
cfg.nSubjects = length(T1_epil_measures_upted.SubjectID);
[counts, pVals, ORvals, hFig] = ecpfunc_plotDiscordantStackedBar(cfg);


cfg = [];
cfg.bestResultsTable = bestResultsTable;
cfg.myCategorical = T1_epil_measures_upted.TLEside;
cfg.discordColumn = 'Gross_Discord_Subs';
cfg.categoryList  = {'Left','Right','Bilateral'};
cfg.doFisher = true;
cfg.binA = {'Bilateral','Right'};
cfg.binB = {'Left'};
cfg.title = 'TLE side';
cfg.nSubjects = length(T1_epil_measures_upted.SubjectID);
[counts, pVals, ORvals, hFig] = ecpfunc_plotDiscordantStackedBar(cfg);


cfg = [];
cfg.bestResultsTable = bestResultsTable;
cfg.myCategorical = T1_epil_measures_upted.AEDcat;
cfg.discordColumn = 'Gross_Discord_Subs';
cfg.categoryList  = {'1','2','3plus'};
cfg.doFisher = true;
cfg.binA = {'1','2'};
cfg.binB = {'3plus'};
cfg.title = 'AED';
cfg.nSubjects = length(T1_epil_measures_upted.SubjectID);
[counts, pVals, ORvals, hFig] = ecpfunc_plotDiscordantStackedBar(cfg);


cfg = [];
cfg.bestResultsTable = bestResultsTable;
cfg.myCategorical = T1_epil_measures_upted.LTGTCcat;
cfg.discordColumn = 'Gross_Discord_Subs';
cfg.categoryList  = {'0','1-5','6-20','21plus'};
cfg.doFisher = true;
cfg.binA = {'0','1-5','6-20'};
cfg.binB = {'21plus'};
cfg.title = 'LTGT';
cfg.nSubjects = length(T1_epil_measures_upted.SubjectID);
[counts, pVals, ORvals, hFig] = ecpfunc_plotDiscordantStackedBar(cfg);

cfg = [];
cfg.bestResultsTable = bestResultsTable;
cfg.myCategorical = T1_epil_measures_upted.SGcat;
cfg.discordColumn = 'Gross_Discord_Subs';
cfg.categoryList  = {'1to2','0','3plus'};
cfg.doFisher = true;
cfg.binA = {'1to2','0'};
cfg.binB = {'3plus'};
cfg.title = 'SG';
cfg.nSubjects = length(T1_epil_measures_upted.SubjectID);
[counts, pVals, ORvals, hFig] = ecpfunc_plotDiscordantStackedBar(cfg);

cfg = [];
cfg.bestResultsTable = bestResultsTable;
cfg.myCategorical = T1_epil_measures_upted.cp_freq_cat;
cfg.discordColumn = 'Gross_Discord_Subs';
cfg.categoryList  = {'1to5','11plus','6to10'};
cfg.doFisher = true;
cfg.binA = {'1to5'};
cfg.binB = {'11plus', '6to10'};
cfg.title = 'cp freq';
cfg.nSubjects = length(T1_epil_measures_upted.SubjectID);
[counts, pVals, ORvals, hFig] = ecpfunc_plotDiscordantStackedBar(cfg);


%%
clc, close all

% --- 1) Plot TLE Side
run_plot_TLEside;  % Creates TLE side stacked bar
doPlotExport(plot_option, save_dir, ...
    sprintf('TLEside_%s_%s', roi, method), 'svg');
disp('--------')

% --- 2) Plot IQ
run_plotIQ;        % Creates IQ stacked bar
doPlotExport(plot_option, save_dir, ...
    sprintf('IQ_%s_%s', roi, method), 'svg');
disp('--------')


% --- 3) Plot EHQ (handedness)
run_plot_EHQ;      % Creates EHQ stacked bar
doPlotExport(plot_option, save_dir, ...
    sprintf('EHQ_%s_%s', roi, method), 'svg');
disp('--------')

% --- 4) Plot AEDcount
run_plot_AED;
% Creates AED count stacked bar
doPlotExport(plot_option, save_dir, ...
    sprintf('AED_%s_%s', roi, method), 'svg');
disp('--------')


run_plot_AED_median
% Creates AED count stacked bar
doPlotExport(plot_option, save_dir, ...
    sprintf('AED_median_%s_%s', roi, method), 'svg');
disp('--------')


% --- 5) Plot LTGTC
run_plot_LTGTC;    % LTGTC stacked bar (multi-level)
doPlotExport(plot_option, save_dir, ...
    sprintf('LTGTC_%s_%s', roi, method), 'svg');
disp('--------')

run_plot_LTGTC_median
doPlotExport(plot_option, save_dir, ...
    sprintf('LTGTC_median_%s_%s', roi, method), 'svg');
disp('--------')

% --- 6) Plot SG frequency
run_plot_SGfreq_median
doPlotExport(plot_option, save_dir, ...
    sprintf('SGfreq_%s_%s', roi, method), 'svg');
disp('--------')

run_plot_SGfreq_median
doPlotExport(plot_option, save_dir, ...
    sprintf('SGfreq_median_%s_%s', roi, method), 'svg');
disp('--------')

% --- 7) Plot CP frequency
run_plot_CPFreq;   % CP_freq_cat (complex partial freq) stacked bar
doPlotExport(plot_option, save_dir, ...
    sprintf('CPfreq_%s_%s', roi, method), 'svg');
disp('--------')

run_plot_CPFreq_median
doPlotExport(plot_option, save_dir, ...
    sprintf('CPfreq_median_%s_%s', roi, method), 'svg');
disp('--------')

%%


lm = fitlm(T, 'Diff_LI ~ TLEside*CP_freq + SG_freq + AEDCount + FSIQ');






