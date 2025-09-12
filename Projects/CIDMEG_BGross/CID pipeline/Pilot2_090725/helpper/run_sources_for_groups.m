function out = run_sources_for_groups(subjectDir, varargin)
% Run group-wise source analysis for folders under a subject:
%   subjectDir/1_slow, subjectDir/2_fast, subjectDir/3_timeout
%
% Steps per group:
%   1) compute noise/data covariance on trials
%   2) average trials
%   3) optional time offset
%   4) compute LCMV/NAI sources on the average
%
% Requirements:
%   - Brainstorm running in MATLAB path (bst_process etc.)
%   - Subject has a default head model selected (GUI)

% -------- options --------
p = inputParser;
addParameter(p,'GroupNames',{'slow','fast','timeout'});
addParameter(p,'FilePrefix','data_2_');         % e.g., 'data_2_'
addParameter(p,'BaseWin',[-0.2 0.0]);           % seconds, for noise covariance
addParameter(p,'DataWin',[0.30 0.90]);          % seconds, for data covariance
addParameter(p,'TimeOffset',0.016);             % seconds; set 0 to skip
addParameter(p,'DoContrast',true);              % Fast - Slow at source level
% Inverse (LCMV / NAI)
inv.Comment        = 'LCMV-NAI (group avg)';
inv.InverseMethod  = 'lcmv';
inv.InverseMeasure = 'nai';
inv.SourceOrient   = {{'fixed'}};
inv.Loose          = 0.2;
inv.UseDepth       = 1;
inv.WeightExp      = 0.5;
inv.WeightLimit    = 10;
inv.NoiseMethod    = 'median';
inv.NoiseReg       = 0.1;
inv.SnrMethod      = 'rms';
inv.SnrRms         = 1e-06;
inv.SnrFixed       = 3;
inv.ComputeKernel  = 1;
inv.DataTypes      = {{'MEG GRAD','MEG MAG'}};
addParameter(p,'Inverse',inv);
parse(p,varargin{:});
opt = p.Results;

% -------- sanity --------
assert(isfolder(subjectDir), 'Subject folder not found: %s', subjectDir);

% -------- run per group --------
bst_report('Start', {});
out = struct();
for gi = 1:numel(opt.GroupNames)
    gName = opt.GroupNames{gi};
    gDir  = fullfile(subjectDir, gName);
    if ~isfolder(gDir)
        fprintf('SKIP group (folder missing): %s\n', gDir);
        continue;
    end

    % Collect trial files in this group
    trialFiles = find_trials_in_group(gDir, opt.FilePrefix);
    if isempty(trialFiles)
        fprintf('SKIP group (no trials): %s\n', gName);
        continue;
    end
    fprintf('Group %s: %d trials\n', gName, numel(trialFiles));

    % 1) Noise covariance (baseline)
    sNoise = bst_process('CallProcess', 'process_noisecov', trialFiles, [], ...
        'baseline',       opt.BaseWin, ...
        'datatimewindow', [], ...
        'sensortypes',    'MEG', ...
        'target',         1, ...  % 1 = noise covariance (baseline)
        'dcoffset',       1, ...
        'identity',       0, ...
        'copycond',       0, ...
        'copysubj',       0, ...
        'copymatch',      0, ...
        'replacefile',    1);

    % 2) Data covariance (active window)
    sData  = bst_process('CallProcess', 'process_noisecov', trialFiles, [], ...
        'baseline',       [], ...
        'datatimewindow', opt.DataWin, ...
        'sensortypes',    'MEG', ...
        'target',         2, ...  % 2 = data covariance (active)
        'dcoffset',       1, ...
        'identity',       0, ...
        'copycond',       0, ...
        'copysubj',       0, ...
        'copymatch',      0, ...
        'replacefile',    1);

    % 3) Average trials
    sAvg = bst_process('CallProcess', 'process_average', trialFiles, [], ...
        'avgtype',         1, ...  % 1 = Everything
        'avg_func',        1, ...  % mean
        'weighted',        0, ...
        'scalenormalized', 0);

    % Optional time offset on the average
    sAvgFinal = sAvg;
    if opt.TimeOffset ~= 0
        sAvgFinal = bst_process('CallProcess', 'process_timeoffset', sAvg, [], ...
            'info',      [], ...
            'offset',    opt.TimeOffset, ...
            'overwrite', 0);
    end

    % 4) Compute sources on the average (uses subject default head model)
    sSrc = bst_process('CallProcess', 'process_inverse_2018', sAvgFinal, [], ...
        'output',     2, ...     % 2 = Full results (one per file)
        'inverse',    opt.Inverse, ...
        'sensortypes','MEG', ...
        'overwrite',  1);

    % Save per-group outputs
    out.(gName).trials   = trialFiles;
    out.(gName).noiseCov = sNoise;
    out.(gName).dataCov  = sData;
    out.(gName).avg      = sAvgFinal;
    out.(gName).sources  = sSrc;
end

% Optional: Fast - Slow contrast at source level (if both exist)
if opt.DoContrast ...
   && isfield(out,'2_fast') && isfield(out,'1_slow') ...
   && ~isempty(out.('2_fast').sources) && ~isempty(out.('1_slow').sources)

    sFast = out.('2_fast').sources;
    sSlow = out.('1_slow').sources;

    % Math: Fast minus Slow
    sDiff = bst_process('CallProcess', 'process_math', sFast, sSlow, ...
        'operator',   'minus', ...
        'timewindow', [], ...
        'isabs',      0, ...
        'comment',    'Fast - Slow (sources)');

    out.contrast.FastMinusSlow = sDiff;
end

bst_report('Save', []);

end

% ----------------- helpers -----------------

function trialFiles = find_trials_in_group(gDir, prefix)
% Accept both data_2_runXX_trialNNN.mat and generic data_2_trialNNN.mat
A = dir(fullfile(gDir, sprintf('%srun*_trial*.mat', prefix)));
B = dir(fullfile(gDir, sprintf('%strial*.mat',     prefix)));
L = [A; B];
L = unique(fullfile({L.folder}, {L.name}), 'stable');
trialFiles = L(:)';
end
