% ICA cleanup using distributed full-rate training windows.
% MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Memory-bounded full-resolution revision: 2026-08-17
%
% ICA weights are estimated from representative windows without temporal
% downsampling. Selected artifact components are subtracted from the full
% recording in blocks. Output sampling rate and duration are unchanged.

%% Setup
clear; clc; close all;
restoredefaultpath;

indir = '/MEG_data/epilepsy';
funcpath = '/MEG_data/MCW_pipeline/Preprocess/func';
mnepath = '/MEG_data/Software/FieldTrip/fieldtrip_2022';

addpath(funcpath);
addpath(fileparts(mfilename('fullpath')), '-begin');

cfgInit = [];
cfgInit.path_tools = '/MEG_data/Software/FieldTrip';
allpath = do_init(cfgInit);

disp('/MEG_data/epilepsy/xx/240612');
cd(indir);
subjdir = uigetdir(indir, 'Choose subject folder');
if isequal(subjdir, 0)
    disp('Operation cancelled.');
    return;
end
cd(subjdir);

%% Options
flag = [];
flag.mnebrowse = 0;
flag.customfilename = 0;

% Full-rate, low-memory ICA settings. Increase training duration only if
% the decomposition is unstable; doing so raises peak memory use.
icaOptions = [];
icaOptions.nWindows = 8;
icaOptions.windowSec = 20;
icaOptions.cleanBlockSec = 20;
icaOptions.maxAutoComponents = 100;
icaOptions.rankTolerance = 1e-8;
icaOptions.rankProbeSamples = 50000;
icaOptions.plotComponents = true;
icaOptions.browseComponents = true;

%% Select data condition
dataConditions = {'raw', '(t) sss', 'all'};
for idx = 1:numel(dataConditions)
    fprintf('%d: %s\n', idx, dataConditions{idx});
end
dcon = input('Select data condition (1=raw, 2=(t)SSS, 3=all): ');
if ~ismember(dcon, 1:3)
    error('Data condition must be 1, 2, or 3.');
end

switch dcon
    case 1
        d = rdir([subjdir, '/**/*.fif']);
        normalizedNames = replace(string({d.name}), '\', '/');
        d = d(~contains(lower(normalizedNames), '/sss/'));
    case 2
        d = rdir([subjdir, '/**/sss/**/*.fif']);
    otherwise
        d = rdir([subjdir, '/**/*.fif']);
end

names = string({d.name});
exclude = endsWith(lower(names), '-eve.fif') | ...
          endsWith(lower(names), '-proj.fif') | ...
          endsWith(lower(names), '_ic_raw.fif') | ...
          endsWith(lower(names), '_ic.fif');
d = d(~exclude);

if isempty(d)
    error('No eligible FIF files were found under %s.', subjdir);
end

confirmed = false;
while ~confirmed
    fprintf('\nAvailable files:\n');
    for i = 1:numel(d)
        fprintf('%4d: %s\n', i, d(i).name);
    end

    fileIndices = input( ...
        'Enter file numbers (e.g., [1 3 5], or 1:length(d)): ');
    if isempty(fileIndices) || any(fileIndices < 1) || ...
            any(fileIndices > numel(d)) || any(mod(fileIndices,1) ~= 0)
        fprintf('Invalid selection. Please try again.\n');
        continue;
    end
    fileIndices = unique(fileIndices, 'stable');

    fprintf('\nSelected files:\n');
    for i = fileIndices
        fprintf('%4d: %s\n', i, d(i).name);
    end

    confirm = input('Confirm files (Yes=1, No=0, Cancel=2): ');
    if isequal(confirm, 1)
        confirmed = true;
    elseif isequal(confirm, 2)
        disp('Operation cancelled.');
        return;
    end
end

requestedComponents = input([ ...
    'Enter ICA component count (0 or [] = estimate effective rank): ']);
if isempty(requestedComponents)
    requestedComponents = 0;
end
if ~isscalar(requestedComponents) || requestedComponents < 0 || ...
        mod(requestedComponents, 1) ~= 0
    error('ICA component count must be zero or a positive integer.');
end

layoutCfg = [];
layoutCfg.layout = 'neuromag306mag.lay';
neuromagLayout = ft_prepare_layout(layoutCfg);

%% Process selected files
for fileNumber = fileIndices
    datafile = d(fileNumber).name;
    fprintf('\n============================================================\n');
    fprintf('Processing: %s\n', datafile);
    fprintf('============================================================\n');

    switch dcon
        case 1
            f_data = do_read_neuromag_raw_nofilter(datafile);
        otherwise
            readCfg = [];
            readCfg.channel = {'MEG'};
            readCfg.datafile = datafile;
            f_data = ft_preprocessing(readCfg);
    end

    fprintf('Loaded %d channels x %d samples at %.3f Hz.\n', ...
        numel(f_data.label), size(f_data.trial{1},2), f_data.fsample);
    originalSampleCount = size(f_data.trial{1},2);
    originalSamplingRate = double(f_data.fsample);

    lowmemCfg = icaOptions;
    lowmemCfg.layout = neuromagLayout;
    lowmemCfg.numcomponent = requestedComponents;

    [cln_data, bic, icaInfo] = ...
        do_ica_fullrate_windowed(lowmemCfg, f_data);

    % Release the input reference before diagnostics/export. cln_data keeps
    % the complete, full-sampling-rate cleaned recording.
    clear f_data

    if size(cln_data.trial{1},2) ~= originalSampleCount || ...
            abs(double(cln_data.fsample)-originalSamplingRate) > ...
            10*eps(originalSamplingRate)
        error('Full-rate cleanup changed the sample count or sampling rate.');
    end
    fprintf('Verified full resolution: %d samples at %.3f Hz.\n', ...
        originalSampleCount, originalSamplingRate);

    fprintf('ICA rank estimate: %d | components learned: %d\n', ...
        icaInfo.rankEstimate, icaInfo.numcomponent);
    fprintf('Training duration: %.1f seconds at %.3f Hz\n', ...
        icaInfo.trainingDurationSec, icaInfo.fsample);

    [rmsMag, rmsGrad] = calculate_meg_rms_blockwise(cln_data, 20);
    fprintf('RMS mag:  %.3g\n', rmsMag);
    fprintf('RMS grad: %.3g\n', rmsGrad);
    fprintf('Ratio grad/mag: %.3g\n', rmsGrad/rmsMag);

    if isempty(bic)
        fprintf('No components selected; no output FIF was written.\n');
        clear cln_data icaInfo
        continue;
    end

    %% Export full-rate cleaned data
    addpath(mnepath);
    addpath(funcpath);

    [savedir, inputBase, inputExt] = fileparts(datafile);
    if endsWith(lower(inputBase), '_raw')
        outputBase = [inputBase(1:end-4), '_ic_raw'];
    else
        outputBase = [inputBase, '_ic'];
    end
    outfile = fullfile(savedir, [outputBase, inputExt]);

    fprintf('Proposed output: %s\n', outfile);
    if flag.customfilename == 1
        nameOkay = input('Name looking okay (Yes=1, No=0)? ');
        if isequal(nameOkay, 0)
            newName = input('Enter new filename including .fif: ', 's');
            outfile = fullfile(savedir, newName);
        end
    end

    if isfile(outfile)
        overwriteAnswer = input( ...
            sprintf('Output exists. Overwrite %s? (Yes=1/No=0): ', outfile));
        if ~isequal(overwriteAnswer, 1)
            fprintf('Skipped existing output.\n');
            clear cln_data icaInfo
            continue;
        end
    end

    exportCfg = [];
    exportCfg.infile = datafile;
    exportCfg.outfile = outfile;
    exportCfg.cln_data = cln_data;

    switch dcon
        case 1
            helperfif = input( ...
                'Enter helper FIF file (full path): ', 's');
            if ~isfile(helperfif)
                error('Helper FIF file not found: %s', helperfif);
            end
            exportCfg.helper = helperfif;
            do_mne_ex_read_write_raw_helper(exportCfg);
        otherwise
            do_mne_ex_read_write_raw_SI_ver2(exportCfg);
    end

    fprintf('Completed. Full-rate cleaned data are ready for review:\n');
    fprintf('  %s\n', outfile);

    if flag.mnebrowse == 1
        system(sprintf('mbrowse "%s" &', outfile));
    end

    clear cln_data icaInfo
end

disp('Script execution completed.');
close all;

function [rmsMag, rmsGrad] = calculate_meg_rms_blockwise(data, blockSec)
mag = ft_chantype(data.label, 'megmag');
grad = ft_chantype(data.label, 'megplanar');
nSamples = size(data.trial{1}, 2);
blockSamples = max(1, round(blockSec*data.fsample));

magSumSquares = 0;
magCount = 0;
gradSumSquares = 0;
gradCount = 0;

for firstSample = 1:blockSamples:nSamples
    lastSample = min(firstSample+blockSamples-1, nSamples);
    index = firstSample:lastSample;

    x = double(data.trial{1}(mag, index));
    magSumSquares = magSumSquares + sum(x(:).^2, 'omitnan');
    magCount = magCount + nnz(isfinite(x));

    x = double(data.trial{1}(grad, index));
    gradSumSquares = gradSumSquares + sum(x(:).^2, 'omitnan');
    gradCount = gradCount + nnz(isfinite(x));
end

rmsMag = sqrt(magSumSquares/max(1,magCount));
rmsGrad = sqrt(gradSumSquares/max(1,gradCount));
end
