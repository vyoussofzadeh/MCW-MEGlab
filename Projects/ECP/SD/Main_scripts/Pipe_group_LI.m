
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
data_save_dir = '/data/MEG/Research/ECP/Semantic_Decision/Results_prestim/';
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
%         cfg.datatag = 'wDICS_hshape';
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

% cfg = [];
% cfg.thre = fMRI_thre; cfg.LI = fmri_LIs_val;
% fmri_LIs_trn = do_ternary_classification2(cfg);
% size(fmri_LIs_trn);

switch LI_analysis
    case 3
        Pt_ID  = LI_pt.Magnitude.setup.S_data{1}.sFiles_subid(IA)';
        Pt_sfiles  = LI_pt.Magnitude.setup.S_data{1}.sFiles_in(IA)';
    otherwise
        Pt_ID  = LI_pt.Magnitude.setup.S_data.sFiles_subid(IA)';
        Pt_sfiles  = LI_pt.Magnitude.setup.S_data.sFiles_in(IA)';
end

%%
fmri_LIs_ROIs = [fmri_LIs.val.language_Angular, fmri_LIs.val.language_Frontal, ...
    fmri_LIs.val.language_Temporal, fmri_LIs.val.language_Lateral];
fmri_LIs_ROIs = fmri_LIs_ROIs(IB_megfmri,:);

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
run_plot_fixedinterval_corcon_rois3

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

%%  fMRI & MEG LI dominance summary
Dominance  =[];
for i = 1:4
    tmp = bestResultsTable.fMRI_LI_tern{i};
    Dominance.(bestResultsTable.ROI{i}).fmri = [length(find(tmp == 1)), length(find(tmp == 0)), length(find(tmp == -1))];
    Dominance.(bestResultsTable.ROI{i}).fmri_perc = Dominance.(bestResultsTable.ROI{i}).fmri./ length(tmp) .* 100;

    tmp = bestResultsTable.MEG_LI_tern{i};
    Dominance.(bestResultsTable.ROI{i}).meg = [length(find(tmp == 1)), length(find(tmp == 0)), length(find(tmp == -1))];
    Dominance.(bestResultsTable.ROI{i}).meg_perc = Dominance.(bestResultsTable.ROI{i}).meg./ length(tmp) .* 100;
    
    Dominance.(bestResultsTable.ROI{i}).method = bestResultsTable.Best_LI_Method{i};
    Dominance.(bestResultsTable.ROI{i}).side = {'Left';'Bilateral';'Right'};

end

%% Task and epil. measures
%- Response (reaction) time
run_responsereaction

% Task Performance
run_taskperformance

%- Epilepsy Metrics
T1_epil_measures = ecpfunc_read_epil_measures();   % Load additional measures

%% Noise SNR
plot_option = 1;
snrDataDir  = '/data/MEG/Research/MEGIN project/Scripts/SNR/SNR_Data';
save_dir_snr    = '/data/MEG/Research/MEGIN project/Scripts/SNR/Plots';

T_snr = ecp_func_noiseSNRAnalysis2(plot_option, snrDataDir, save_dir_snr);
disp(T_snr);

%% Combine all variables
combined1 = outerjoin(T_patn_MEGfMRI, accuracyResults_updt, 'Keys', 'SubjectID', 'MergeKeys', true);

% Merge the result with T1_epil_measures_tbl
final_combined = outerjoin(combined1, T1_epil_measures, 'LeftKeys', 'SubjectID', 'RightKeys', 'SubjectID', 'MergeKeys', true);
[commonSubs, idxFinal, idxPt] = intersect(final_combined.SubjectID, Pt_ID);
final_combined_updt       = final_combined(idxFinal,:);

%
% Make sure SubjectID columns match in format and case
if iscell(final_combined_updt.SubjectID)
    final_combined_updt.SubjectID = upper(string(final_combined_updt.SubjectID));
else
    final_combined_updt.SubjectID = upper(final_combined_updt.SubjectID);
end
if iscell(T_snr.SubjectID)
    T_snr.SubjectID = upper(string(T_snr.SubjectID));
else
    T_snr.SubjectID = upper(T_snr.SubjectID);
end

% Perform a LEFT JOIN using outerjoin with 'Type','left'
%    - This keeps ALL subjects from final_combined_updt
%    - Adds columns from T_snr where SubjectID matches
%    - If a subject is missing in T_snr, those columns become NaN/empty
final_combined_snr = outerjoin(final_combined_updt, T_snr, ...
    'LeftKeys','SubjectID', 'RightKeys','SubjectID', ...
    'Type','left', ...           % left join
    'MergeKeys',true);

% Inspect the result
head(final_combined_snr)

%%  fMRI & MEG LI dominance summary
Dominance  =[];
for i = 1:4
    tmp = bestResultsTable.fMRI_LI_tern{i};
    Dominance.(bestResultsTable.ROI{i}).fmri = [length(find(tmp == 1)), length(find(tmp == 0)), length(find(tmp == -1))];
    Dominance.(bestResultsTable.ROI{i}).fmri_perc = Dominance.(bestResultsTable.ROI{i}).fmri./ length(tmp) .* 100;

    tmp = bestResultsTable.MEG_LI_tern{i};
    Dominance.(bestResultsTable.ROI{i}).meg = [length(find(tmp == 1)), length(find(tmp == 0)), length(find(tmp == -1))];
    Dominance.(bestResultsTable.ROI{i}).meg_perc = Dominance.(bestResultsTable.ROI{i}).meg./ length(tmp) .* 100;
   
    Dominance.(bestResultsTable.ROI{i}).method = bestResultsTable.Best_LI_Method{i};
    Dominance.(bestResultsTable.ROI{i}).side = {'Left';'Bilateral';'Right'};

end

%% Plot noise SNR
clc
close all

plot_option = 1;
save_dir_test = '/path/to/export';

% ecp_plot_noiseSNR(final_combined_snr, plot_option, save_dir);

ecp_plot_noiseSNRCombined(final_combined_snr, plot_option, save_dir_test)

% plot_option = 1;
% run_plotSNR

%% Investigate Discordant Samples of Best Results and Obtain Corresponding MEG_LI and fMRI_LI
clc
close all

for roi = 3:3
    cfg = [];
    cfg.roi_sel = roi; %lateral, n=3
    cfg.wi = wi;
    cfg.bounds = bounds;
    cfg.plot_option  = 1;
    cfg.save_dir = save_dir;
    cfg.final_combined_updt = final_combined_snr;
    cfg.T1_epil_measures = T1_epil_measures;
    switch opt_sel
        case 1
            cfg.bestResultsTable = bestResultsTable;
        case 4
            cfg.bestResultsTable = bestResultsTable_fixed;
    end
    ecpfunc_assess_gross_discondances2(cfg)
end

%% Beta
T_snr.tSSS_3cat_beta = defineTertiles(T_snr.nSNR_Beta_MEGnetvstSSS);
T_snr.megnet_3cat_beta = defineTertiles(T_snr.nSNR_Beta_tSSSvsRaw);
head(T_snr(:, ["nSNR_Beta_MEGnetvstSSS","nSNR_Beta_tSSSvsRaw"]))

T_snr.tSSS_2cat_beta = defineTwoCategories(T_snr.nSNR_Beta_MEGnetvstSSS);
T_snr.megnet_2cat_beta = defineTwoCategories(T_snr.nSNR_Beta_tSSSvsRaw);

%% Broadband
T_snr.tSSS_3cat_broad = defineTertiles(T_snr.nSNR_Broad_MEGnetvstSSS);
T_snr.megnet_3cat_broad = defineTertiles(T_snr.nSNR_Broad_tSSSvsRaw);
head(T_snr(:, ["nSNR_Broad_MEGnetvstSSS","nSNR_Broad_tSSSvsRaw"]))

T_snr.tSSS_2cat_broad = defineTwoCategories(T_snr.nSNR_Broad_MEGnetvstSSS);
T_snr.megnet_2cat_broad = defineTwoCategories(T_snr.nSNR_Broad_tSSSvsRaw);

%% Fisher analysis
%- Epilepsy Metrics
% Convert or define categorical variables
T1_epil_measures_upted = [];
T1_epil_measures_upted.EHQcat       = defineHandedness(final_combined_updt.EHQ, 40);
T1_epil_measures_upted.IQcat        = defineIQbins(final_combined_updt.NP1WASI_FSIQ);
T1_epil_measures_upted.TLEside      = defineTLEside(final_combined_updt.TLEside);
T1_epil_measures_upted.AEDcat       = defineAEDbins(final_combined_updt.AEDCount);
T1_epil_measures_upted.LTGTCcat     = defineLTGTCbins(final_combined_updt.LTGTC);

T1_epil_measures_upted.SGcat        = defineSGfreqbins(final_combined_updt.SG_freq);

T1_epil_measures_upted.cp_freq_cat  = defineCPfreqbins(final_combined_updt.CP_freq);

T1_epil_measures_upted.AnimalRTcat   = defineRTbins(final_combined_updt.Animal_RT);
T1_epil_measures_upted.SymbolRTcat   = defineRTbins(final_combined_updt.Symbol_RT);
T1_epil_measures_upted.AnimalACCcat  = defineACCbins(final_combined_updt.Animal_ACC);
T1_epil_measures_upted.SymbolACCcat  = defineACCbins(final_combined_updt.Symbol_ACC);
T1_epil_measures_upted.SubjectID     = final_combined_updt.SubjectID;

%% Braodband
% 1) Find matching row indices
[commonSubj, idxInT1, idxInSNR] = intersect( ...
        T1_epil_measures_upted.SubjectID, ...
        T_snr.SubjectID, 'stable');

% 2) Copy the new columns for those matching subjects
T1_epil_measures_upted.tSSS_3cat_broad(idxInT1)   = T_snr.tSSS_3cat_broad(idxInSNR);
T1_epil_measures_upted.megnet_3cat_broad(idxInT1) = T_snr.megnet_3cat_broad(idxInSNR);
T1_epil_measures_upted.tSSS_2cat_broad(idxInT1)   = T_snr.tSSS_2cat_broad(idxInSNR);
T1_epil_measures_upted.megnet_2cat_broad(idxInT1) = T_snr.megnet_2cat_broad(idxInSNR);

T1_epil_measures_upted.tSSS_3cat_beta(idxInT1)   = T_snr.tSSS_3cat_beta(idxInSNR);
T1_epil_measures_upted.megnet_3cat_beta(idxInT1) = T_snr.megnet_3cat_beta(idxInSNR);
T1_epil_measures_upted.tSSS_2cat_beta(idxInT1)   = T_snr.tSSS_2cat_beta(idxInSNR);
T1_epil_measures_upted.megnet_2cat_beta(idxInT1) = T_snr.megnet_2cat_beta(idxInSNR);

%%
clc, close all
% run_fisheranalysis_2x2_beta
% run_fisheranalysis_2x2_broad
run_fisheranalysis_2x2_all_2

% run_fisheranalysis_2x3_beta
% run_fisheranalysis_2x3_broad
run_fisheranalysis_2x3_all_2
cd(save_dir)



