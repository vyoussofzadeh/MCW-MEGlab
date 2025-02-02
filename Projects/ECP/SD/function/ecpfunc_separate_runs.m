function S_data_sep = ecpfunc_separate_runs(S_data_in)
% ECPFUNC_SEPARATE_RUNS
%   Splits the sFiles_2 and sFiles_3 fields of S_data_in into Run 1 and Run 2
%   based on the "_runX_" portion of each file name. Also extracts subject
%   IDs if needed.
%
%   INPUT:
%       S_data_in    - a struct with fields:
%                      .sFiles_2 : cell array of file paths (condition 2)
%                      .sFiles_3 : cell array of file paths (condition 3)
%                      .subjs_2  : cell array of subject IDs for sFiles_2
%                      .subjs_3  : cell array of subject IDs for sFiles_3
%
%   OUTPUT:
%       S_data_sep   - a struct with separate fields for Run1 and Run2:
%                      .sFiles_2_run1, .sFiles_2_run2, etc.
%                      and similarly subjs_2_run1, subjs_2_run2, etc.
%
%   Example:
%       S_data_sep = ecpfunc_separate_runs(S_data);

    % Initialize the output struct
    S_data_sep = struct();

    % --- Separate Condition 2 (sFiles_2) ---
    [S_data_sep.sFiles_2_run1, ...
     S_data_sep.sFiles_2_run2, ...
     S_data_sep.subjs_2_run1, ...
     S_data_sep.subjs_2_run2] = separate_by_run(S_data_in.sFiles_2, S_data_in.subjs_2);

    % --- Separate Condition 3 (sFiles_3) ---
    [S_data_sep.sFiles_3_run1, ...
     S_data_sep.sFiles_3_run2, ...
     S_data_sep.subjs_3_run1, ...
     S_data_sep.subjs_3_run2] = separate_by_run(S_data_in.sFiles_3, S_data_in.subjs_3);

end

% -------------------------------------------------------------------------
% Helper function: separate_by_run
% -------------------------------------------------------------------------
function [files_run1, files_run2, subs_run1, subs_run2] = separate_by_run(fileList, subList)
% Splits a single fileList + subList into Run1 and Run2 based on a "_runX_" pattern.

    % Initialize outputs
    files_run1 = {};
    files_run2 = {};
    subs_run1  = {};
    subs_run2  = {};

    for i = 1:length(fileList)
        fPath = fileList{i};
        runNum = detect_run(fPath);

        if runNum == 1
            files_run1{end+1} = fPath;
            subs_run1{end+1}  = subList{i};
        elseif runNum == 2
            files_run2{end+1} = fPath;
            subs_run2{end+1}  = subList{i};
        else
            % If you need to handle more than 2 runs,
            % you could adjust or add more conditions here.
            warning('File %s has runNum=%d, ignoring in this example.', fPath, runNum);
        end
    end
end

% -------------------------------------------------------------------------
% Helper function: detect_run
% -------------------------------------------------------------------------
function runNum = detect_run(filePath)
% Detect run number using a regex on the filename, e.g., "_run1_" or "_run2_"
    runPattern = '_run(\d+)_'; % looks for e.g. "_run1_" or "_run2_"
    tokens = regexp(filePath, runPattern, 'tokens');
    if ~isempty(tokens)
        runNum = str2double(tokens{1}{1});
    else
        % Default if not found
        runNum = 1; 
    end
end
