%% The CID (Pilot2) MEG

% Performance (acc) analysis 1
% Writtern by MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Lastest update: 09/12/2025

cd /group/bgross/work/CIDMEG/Analysis/BrainstormProcess/Brainstorm_db/Pilot2_092025/data/mcwa086_v1/mcwa086_v1_Run_02
addpath('/group/bgross/work/CIDMEG/Analysis/Pilot2_090725/')

% Average RT=fast trials in Run 2
avg_fast = average_by_tag(2, 'RT=fast', 'data_2_avg_RTfast.mat');

% Average RT=slow trials in Run 2
avg_slow = average_by_tag(2, 'RT=slow', 'data_2_avg_RTslow.mat');

% (Optional) Average Key=1 trials in Run 2
avg_key1 = average_by_tag(2, 'Key=1', 'data_2_avg_Key1.mat');


% You already have T in workspace (from your CSV/table)
subjectDir = '/group/bgross/work/CIDMEG/Analysis/BrainstormProcess/Brainstorm_db/Pilot2_092025/data/mcwa086_v1';

% Average across ALL runs found in T:
results = avg_by_rt_from_table(T, subjectDir, 'IncludeTimeout', false, 'DoContrast', true);

% T is your behavioral table already in workspace
subjectDir = '/group/bgross/work/CIDMEG/Analysis/BrainstormProcess/Brainstorm_db/Pilot2_092025/data/mcwa086_v1';

run_behav

out = avg_allruns_by_rt_from_table(T, subjectDir, ...
       'IncludeTimeout', false, ...     % exclude timeouts
       'DoContrast', true, ...
       'FilePrefix', 'data_2_');        % you said files start with data_2_
%%
subjectDir = '/group/bgross/work/CIDMEG/Analysis/BrainstormProcess/Brainstorm_db/Pilot2_092025/data/mcwa086_v1';

out = collect_trials_to_groups(T, subjectDir, 'FilePrefix','data_2_');

%%

% Add a text group column that the function can use
T.grp = strings(height(T),1);
T.grp(T.rt_category==1) = "fast";
T.grp(T.rt_category==2) = "slow";
T.grp(T.rt_category==3) = "timeout";

out = collect_trials_to_groups_acc(T, subjectDir, ...
    'FilePrefix','data_2_', ...          % matches your trial files
    'GroupVar','group', ...              % or auto-detect
    'AccVar','is_correct', ...           % logical/0-1 or text
    'FileVar','FileName', ...            % if T already has filenames
    'Action','copy');                    % 'copy' | 'link' | 'none'


%% run source analysis manually (BS GUI) ... having issues ....!!!

% Make sure the subject default head model is set (e.g., from Run_02) in the GUI.
% out = run_sources_for_groups(subjectDir, ...
%     'FilePrefix','data_2_', ...
%     'BaseWin',[-0.2 0.0], ...
%     'DataWin',[0.30 0.90], ...
%     'TimeOffset',0.016, ...
%     'DoContrast',true);
