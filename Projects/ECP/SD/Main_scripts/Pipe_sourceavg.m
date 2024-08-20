%% ECP Semantic decision task dataset, Medical College of Wisconsin

% Script: BS Process (Laterality analysis)
% Project: ECP_SD
% Writtern by: Vahab Youssof Zadeh
% Update: 06/27/2023

clear; clc, close('all'); warning off,

%% Paths
restoredefaultpath
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/run')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/function')
Run_setpath

%%
flags = [];
flags.plot_atlasnetworks = 0;

%% HCP Atlas
src_fname = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/cortex_pial_low.fs';
data_save_dir = '/data/MEG/Research/ECP/Semantic_Decision/Results_prestim';
glass_atlas = '/data/MEG/Vahab/Github/MCW_MEGlab/tools/Atlas/HCP/HCP atlas for Brainstorm/Best/scout_mmp_in_mni_symmetrical_final_updated.mat';
glass_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/Glasser';

cfg = struct('src_fname', src_fname, 'glass_dir', glass_dir, 'glass_atlas', glass_atlas, 'plotflag', 0);
Data_hcp_atlas = ecpfunc_hcp_atlas2(cfg);

%%
% Lateral network
% close all
cfg = [];
cfg.src_fname = src_fname;
cfg.network_sel = [11];
cfg.Data_hcp_atlas = Data_hcp_atlas;
cfg.plotflag = 0;
cfg.fixedcolor = [0,0.7,0];
[idx_L, idx_R, src]  = do_plot_hcp_network(cfg);
net_rois = 'ftp';
disp(Data_hcp_atlas.groups_labels)

%
cfg = []; cfg.idx_L = idx_L; cfg.idx_R = idx_R; cfg.Data_hcp_atlas = Data_hcp_atlas;
cfg.src_fname = src_fname;

% cfg.export = 1; cfg.savedir = fullfile(outdir,'group');
% cfg.network_sel = [1,2,6]; do_map_HCP_net_sel(cfg);
% net_label = 'Fronto_tempro_pri';

save_dir_atlas = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/results/atlas_roi';

colorcode = {[0,114,189]; [217,83,25]; [237,177,32];[126,47,142]; [119,172,48]};
network_sel = [1,2,5,6,11]; % LI networks compared between MEG and fMRI

if flags.plot_atlasnetworks == 1
    
    % Inspecting atlas networks common between MEG and fMRI
    for i = 1:length(colorcode)
        cfg.network_sel = network_sel(i);
        cfg.fixedcolor = colorcode{i}/256;
        do_map_HCP_net_sel(cfg);
        title(Data_hcp_atlas.groups_labels{cfg.network_sel});
        
        cfg2 = [];
        cfg2.outdir = save_dir_atlas;
        filename = Data_hcp_atlas.groups_labels{cfg.network_sel};
        cfg2.filename = filename;
        cfg2.type = 'fig';
        do_export_fig(cfg2)
    end
    
    
    % Inspecting all atlas networks areas (MEG only)
    for i = 1:11
        cfg.network_sel = i;
        cfg.fixedcolor = colorcode/256;
        do_map_HCP_net_sel(cfg);
        title(Data_hcp_atlas.groups_labels{cfg.network_sel});
        
        cfg2 = [];
        cfg2.outdir = save_dir_atlas;
        filename = Data_hcp_atlas.groups_labels{cfg.network_sel};
        cfg2.filename = filename;
        cfg2.type = 'fig';
        do_export_fig(cfg2)
    end
end

%% Run Brainstorm
Run_BS
brainstorm

%%
Run_load_surface_template

%%
LI_analysis_label = {'DICS_indirect','DICS_directcontrast','LCMV_anim_vs_Symb','-','DICS_anim', 'DICS_contrast_prestim', 'dSPM_contrast', ...
    'DICS_directcontrast_onlyprestim'};

for i = 1:length(LI_analysis_label)
    disp([num2str(i) ') ' LI_analysis_label{i}]);
end

LI_analysis = input('');

%%
datadir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data';
BS_data_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/data_full';

cfg = [];
cfg.protocol = protocol;
cfg.datadir = datadir;
cfg.BS_data_dir = BS_data_dir;

switch LI_analysis
    case {1,5}
        cfg.BS_data_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/data';
        %         cfg.datamask = 'wDICS_baseline_18_4_50ms';
        cfg.datatag = 'wDICS_baseline_18_4';
        S_data = ecpfunc_read_sourcemaps_dics(cfg);
    case 2
        %         cfg.datamask = 'wDICS_18_4_50ms';
        %         S_data = ecpfunc_read_sourcemaps_dics_contrast(cfg);
        
        cfg.BS_data_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/data';
        cfg.datatag = 'wDICS_contrast_18_4';
        S_data = ecpfunc_read_sourcemaps_dics_contrast(cfg);
        
    case 3
        cfg.datamask = fullfile('./Group_analysis/LCMV/results_average*.mat');
        S_data = ecpfunc_read_sourcemaps(cfg);
    case 6
        cfg.datamask = 'wDICS_18_4_50ms';
        S_data = ecpfunc_read_sourcemaps_prestim(cfg);
    case 7
        cfg.datamask = 'dSPM_contrast';
        S_data = ecpfunc_read_sourcemaps_contrast(cfg);
    case 8
        cfg.datamask = 'wDICS_18_4_presim500ms';
        S_data = ecpfunc_read_sourcemaps_prestim(cfg);
end












datadir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data';
BS_data_dir = '/data/MEG/Research/ECP/Semantic_Decision/BS_database/data_full';

cfg = [];
cfg.protocol = protocol;
cfg.datadir = datadir;
cfg.BS_data_dir = BS_data_dir;

switch LI_analysis
    case {1,5}
        cfg.datatag = 'wDICS_baseline_18_4_50ms';
        S_data = ecpfunc_read_sourcemaps_dics(cfg);
    case 2
        cfg.datatag = 'wDICS_18_4_50ms';
        S_data = ecpfunc_read_sourcemaps_dics_contrast(cfg);
    case 3
        cfg.datamask = fullfile('./Group_analysis/LCMV/results_average*.mat');
        S_data = ecpfunc_read_sourcemaps(cfg);
%     case 4
%         cfg.datamask = fullfile('./Group_analysis/LCMV/results_abs*.mat');
%         S_data = ecpfunc_read_sourcemaps_contrast(cfg);
    case 6
        cfg.datatag = 'wDICS_contrast_18_4_50ms';
        S_data = ecpfunc_read_sourcemaps_dics_contrast(cfg);
    case 7
        cfg.datatag = 'dSPM_contrast';
        S_data = ecpfunc_read_sourcemaps_contrast(cfg);
end

%% TLE side (PT only)
patn_neuropsych_data = ecpfunc_read_patn_neuropsych();

%% Subject demog details
switch LI_analysis
    case {1,3 ,5}
        cfg = []; cfg.subjs_3 = S_data.subjs_3; cfg.subjs_2 = S_data.subjs_2;
        cfg.sFiles_3 = S_data.sFiles_3; cfg.sFiles_2 = S_data.sFiles_2;
        sub_demog_data = ecpfunc_read_sub_demog(cfg);
    case {2, 6, 7, 8}
        cfg = []; cfg.subjs = S_data.subjs;
        cfg.sFiles = S_data.sFiles_32;
        sub_demog_data = ecpfunc_read_sub_demog_contrast(cfg);
end

%%
switch LI_analysis
    case {1,3}
        
        % Process: Difference: A-B, abs
        sFiles = bst_process('CallProcess', 'process_diff_ab', sub_demog_data.sFiles_anim_patn, sub_demog_data.sFiles_symb_patn, ...
            'source_abs', 1);
        
        % Process: Average: Everything
        bst_process('CallProcess', 'process_average', sFiles, [], ...
            'avgtype',         1, ...  % Everything
            'avg_func',        1, ...  % Arithmetic average:  mean(x)
            'weighted',        0, ...
            'scalenormalized', 0);
        
        
        % Process: Perm t-test paired [0.000s,2.100s]          H0:(A=B), H1:(A<>B)
        bst_process('CallProcess', 'process_test_permutation2p', sub_demog_data.sFiles_anim_patn, sub_demog_data.sFiles_symb_patn, ...
            'timewindow',     [7.21644966e-16, 2.1], ...
            'scoutsel',       {}, ...
            'scoutfunc',      1, ...  % Mean
            'isnorm',         0, ...
            'avgtime',        0, ...
            'iszerobad',      1, ...
            'Comment',        '', ...
            'test_type',      'ttest_paired', ...  % Paired Student's t-test T = mean(A-B) / std(A-B) * sqrt(n)
            'randomizations', 1000, ...
            'tail',           'two');  % Two-tailed
        
        
        timeWindowStart = 0; % Start time in ms
        timeWindowEnd = 1.99; % End time in ms
        interval = .3; % Interval in ms
        
        % Convert ms to seconds for Brainstorm functions
        for startTime = timeWindowStart: interval: timeWindowEnd
            endTime = startTime + interval; % Define the end time of the current window
            % Ensure the time is in seconds for the Brainstorm process
            startTimeSec = startTime;
            endTimeSec = endTime;
            
            disp([startTimeSec, endTimeSec])
            %             pause,
            output = calculateIntervalAverage(brainstormData, startTimeSec, endTimeSec);
            key = ['brainstormData_', [num2str(1000*startTimeSec), '_', num2str(1000*endTimeSec)]];
            eval([key ' = output;']); % Assign the value to the dynamically named variable
            
            pause(2)
        end
    case 5
        
        bst_process('CallProcess', 'process_average', sub_demog_data.sFiles_anim_patn, [], ...
            'avgtype',         1, ...  % Everything
            'avg_func',        1, ...  % Arithmetic average:  mean(x)
            'weighted',        0, ...
            'Comment', 'avg_anim_patn', ...
            'scalenormalized', 0);
        
        % Process: t-test zero [-500ms,1950ms]          H0:(X=0), H1:(X<>0)
        bst_process('CallProcess', 'process_test_parametric1', sub_demog_data.sFiles_anim_patn, [], ...
            'timewindow',    [-0.5, 1.95], ...
            'scoutsel',      {}, ...
            'scoutfunc',     1, ...  % Mean
            'isnorm',        0, ...
            'avgtime',       0, ...
            'Comment',       '', ...
            'Comment', 'stats_anim_patn', ...
            'test_type',     'ttest_onesample', ...  % One-sample Student's t-test    X~N(m,s)t = mean(X) ./ std(X) .* sqrt(n)      df=n-1
            'tail',          'two');  % Two-tailed
        
        
    case {2, 6, 7, 8}
        
%         % Process: Extract time: [0.000s,2.100s]
%         sFiles = bst_process('CallProcess', 'process_extract_time', sub_demog_data.sFiles_patn(1:2), [], ...
%             'timewindow', [7.21644966e-16, 0.3], ...
%             'overwrite',  0);
        
        % PT
        % Process: Average: Everything
        sFiles_mean = bst_process('CallProcess', 'process_average', sub_demog_data.sFiles_patn, [], ...
            'avgtype',         1, ...  % Everything
            'avg_func',        1, ...  % Arithmetic average:  mean(x)
            'weighted',        0, ...
            'Comment', 'avg_patn', ...
            'scalenormalized', 0);
        
%         % Process: t-test zero [-500ms,1950ms]          H0:(X=0), H1:(X<>0)
%         bst_process('CallProcess', 'process_test_parametric1', sub_demog_data.sFiles_patn, [], ...
%             'timewindow',    [-0.5, 1.95], ...
%             'scoutsel',      {}, ...
%             'scoutfunc',     1, ...  % Mean
%             'isnorm',        0, ...
%             'avgtime',       0, ...
%             'Comment',       '', ...
%             'Comment', 'stats_patn', ...
%             'test_type',     'ttest_onesample', ...  % One-sample Student's t-test    X~N(m,s)t = mean(X) ./ std(X) .* sqrt(n)      df=n-1
%             'tail',          'two');  % Two-tailed
        
%         % HC
%         % Process: Average: Everything
%         bst_process('CallProcess', 'process_average', sub_demog_data.sFiles_ctrl, [], ...
%             'avgtype',         1, ...  % Everything
%             'avg_func',        1, ...  % Arithmetic average:  mean(x)
%             'weighted',        0, ...
%             'Comment', 'avg_ctrl', ...
%             'scalenormalized', 0);
%         
%         % Process: t-test zero [-500ms,1950ms]          H0:(X=0), H1:(X<>0)
%         bst_process('CallProcess', 'process_test_parametric1', sub_demog_data.sFiles_ctrl, [], ...
%             'timewindow',    [-0.5, 1.95], ...
%             'scoutsel',      {}, ...
%             'scoutfunc',     1, ...  % Mean
%             'isnorm',        0, ...
%             'avgtime',       0, ...
%             'Comment',       '', ...
%             'Comment', 'stats_ctrl', ...
%             'test_type',     'ttest_onesample', ...  % One-sample Student's t-test    X~N(m,s)t = mean(X) ./ std(X) .* sqrt(n)      df=n-1
%             'tail',          'two');  % Two-tailed
        
end

%%
cfg.strt = -0.3;
cfg.spt = 2;
cfg.overlap = 0;
% cfg.overlap = 0.01;
% cfg.overlap = 0.15;
cfg.linterval = 0.15;
% cfg.linterval = 0.1;
wi  = do_time_intervals(cfg);

sFiles = sFiles_mean;

for startTime = 1:length(wi)
%     endTime = startTime + interval; % Define the end time of the current window
    % Ensure the time is in seconds for the Brainstorm process
    startTimeSec = wi(startTime,1);
    endTimeSec = wi(startTime,2);
    
    disp([startTimeSec, endTimeSec])
    % Process: MEAN: [startTimeSec, endTimeSec], abs
    bst_process('CallProcess', 'process_average_time', sFiles, [], ...
        'timewindow', [startTimeSec, endTimeSec], ...
        'avg_func', 'mean', ...  % Arithmetic average: mean(x)
        'overwrite', 0, ...
        'source_abs', 1);
    
    pause(2)
    
    % Optionally, you can add commands here to save or process the output
end


%%
% % Input files
% sFiles = {'Group_analysis/wDICS_contrast_18_4/results_average_240216_1127.mat'};
% 
% % Define the time window end points in milliseconds
% timeWindowStart = 0; % Start time in ms
% timeWindowEnd = 1.99; % End time in ms
% interval = .3; % Interval in ms
% 
% % sFiles = {'Group_analysis/wDICS_18_4_50ms/results_average_240228_1634.mat'};
% 
% sFiles = {'Group_analysis/wDICS_baseline_18_4/results_average_240216_1402.mat'};
% 
% % Define the time window end points in milliseconds
% timeWindowStart = 0; % Start time in ms
% timeWindowEnd = 1.99; % End time in ms
% interval = .3; % Interval in ms
% 
% sFiles = {'Group_analysis/wDICS_18_4_presim500ms/results_average_240304_1523.mat'};

% sFiles = sFiles_mean;
% 
% timeWindowStart = -0.5; % Start time in ms
% timeWindowEnd = 0; % End time in ms
% interval = .25; % Interval in ms
% 
% 
% 
% %%
% % Convert ms to seconds for Brainstorm functions
% for startTime = timeWindowStart: interval: timeWindowEnd
%     endTime = startTime + interval; % Define the end time of the current window
%     % Ensure the time is in seconds for the Brainstorm process
%     startTimeSec = startTime;
%     endTimeSec = endTime;
%     
%     disp([startTimeSec, endTimeSec])
%     % Process: MEAN: [startTimeSec, endTimeSec], abs
%     bst_process('CallProcess', 'process_average_time', sFiles, [], ...
%         'timewindow', [startTimeSec, endTimeSec], ...
%         'avg_func', 'mean', ...  % Arithmetic average: mean(x)
%         'overwrite', 0, ...
%         'source_abs', 1);
%     
%     pause(2)
%     
%     % Optionally, you can add commands here to save or process the output
% end

% After the loop, you may want to finalize the report or save the last state
% of sFiles, depending on your requirements.
%%
% Define labels for data groups
% Dlabel = {'Groupmean_interval_HC', 'Groupmean_interval_PT'};
%
% % Initialize container for processed data
% newbrainstormData_all = cell(length(Dlabel), 1);
%
% for select_data = 1:length(Dlabel)
%
%     % Configuration for data selection
%     cfg = [];
%     cfg.sub_demog_data = sub_demog_data;
%     cfg.select_data = select_data;
%     S_data_sel = ecpfunc_select_data_contrast(cfg);
%
%     % Preallocate for averaged data
%     tmp_avg = [];
%
%     % Load and average data for each file
%     for i = 1:length(S_data_sel.sFiles_in)
%         tmp = load(fullfile(BS_data_dir, S_data_sel.sFiles_in{i}));
%
%         for j = 1:size(wi,1)
%             timind1 = nearest(tmp.Time, wi(j,1));
%             timind2 = nearest(tmp.Time, wi(j,2));
%             tmp_avg(i,j,:) = nanmean(tmp.ImageGridAmp(:,timind1:timind2), 2);
%         end
%     end
%
%     % Compute the average across all files
%     dics = squeeze(mean(tmp_avg, 1));
%
%     % Prepare the time vector for the averaged data
%     DataMat.Time = wi(1,1):0.01:wi(end,2);
%     numVertices = size(tmp.ImageGridAmp, 1);
%
%     % Initialize the matrix for averaged source data
%     avgSourceData = zeros(numVertices, length(DataMat.Time));
%
%     % Populate the matrix with computed averages
%     for i = 1:size(wi,1)
%         start_idx = find(abs(DataMat.Time - wi(i,1)) < 1e-10, 1);
%         end_idx = find(abs(DataMat.Time - wi(i,2)) < 1e-10, 1);
%         avgSourceData(:, start_idx:end_idx) = repmat(dics(i,:)', 1, end_idx-start_idx+1);
%     end
%
%     % Create a new Brainstorm data structure with the averaged data
%     newbrainstormData = brainstormData;
%     newbrainstormData.ImageGridAmp = avgSourceData;
%     newbrainstormData.Comment = Dlabel{select_data};
%
%     % Store the new data structure in the array
%     newbrainstormData_all{select_data} = newbrainstormData;
% end
%
% % Separate variables for each group for easier access
% newbrainstormData_HC = newbrainstormData_all{1};
% newbrainstormData_PT = newbrainstormData_all{2};




