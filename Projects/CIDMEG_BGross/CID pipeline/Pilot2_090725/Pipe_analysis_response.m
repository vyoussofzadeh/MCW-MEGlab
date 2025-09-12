%% The CID (Pilot2) MEG

% Performance analysis 2 (response 1 vs. 2 contrast)
% Writtern by MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Lastest update: 09/12/2025

addpath('/group/bgross/work/CIDMEG/Analysis/Pilot2_090725/helpper')

% Folder of per-run CSVs like: cidmeg004_times_run1.csv, cidmeg004_times_run2.csv, ...
csvDir = '/group/bgross/work/CIDMEG/Analysis/Pilot2_090725/Event_times';
outDir = fullfile(csvDir, 'resp_split_events');  % will be created

generate_resp_trial_events(csvDir, outDir, 'TimeUnit','s', 'MaxRespLag',10);

%%

T = build_trial_response_table(csvDir, fullfile(csvDir,'trial_response_table.csv'), ...
      'TimeUnit','s','MaxRespLag',10);
%%
% Execute:
bst_bucket_trials_by_response_copy(T, 'mcwa086_v1', 1, ...
    'SourceSuffix','_resp', 'TargetResp1','Resp1', 'TargetResp2','Resp2', ...
    'DryRun', false);

% Execute (copies files + tags comments with [Run XX]):
bst_bucket_trials_by_response_copy(T, 'mcwa086_v1', 1, ...
  'SourceSuffix','_resp', 'TargetResp1','Resp1', 'TargetResp2','Resp2', ...
  'DryRun', false, 'TagRunInComment', true);