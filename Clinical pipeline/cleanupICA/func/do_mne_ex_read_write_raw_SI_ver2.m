function do_mne_ex_read_write_raw_SI_ver2(cfg)
% Stream cleaned MEG data to FIF without loading the full raw file twice.

arguments
    cfg struct
end

requiredFields = {'infile','outfile','cln_data'};
for k = 1:numel(requiredFields)
    if ~isfield(cfg, requiredFields{k}) || isempty(cfg.(requiredFields{k}))
        error('do_mne_ex_read_write_raw_SI_ver2:MissingInput', ...
            'cfg.%s is required.', requiredFields{k});
    end
end
if ~isfield(cfg, 'blockSec') || isempty(cfg.blockSec)
    cfg.blockSec = 20;
end

infile = cfg.infile;
outfile = cfg.outfile;
clnData = cfg.cln_data;

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

raw = fiff_setup_read_raw(infile);
inputCleanup = onCleanup(@() close_fif_safely(raw.fid));
picks = 1:numel(raw.info.chs);

isMeg = false(size(picks));
for k = 1:numel(picks)
    isMeg(k) = raw.info.chs(k).kind == FIFF.FIFFV_MEG_CH;
end
megIndex = find(isMeg);
if isempty(megIndex)
    error('do_mne_ex_read_write_raw_SI_ver2:NoMEGChannels', ...
        'No MEG channels were found in %s.', infile);
end

hdr = ft_read_header(infile);
if numel(hdr.label) ~= numel(picks)
    error('do_mne_ex_read_write_raw_SI_ver2:HeaderMismatch', ...
        'FieldTrip and MNE report different channel counts for %s.', infile);
end
inputMegLabels = hdr.label(isMeg);
[inputSelection, cleanSelection] = ...
    match_str(inputMegLabels, clnData.label);

if numel(inputSelection) ~= numel(inputMegLabels)
    missing = setdiff(inputMegLabels, clnData.label, 'stable');
    error('do_mne_ex_read_write_raw_SI_ver2:MissingCleanChannels', ...
        'Cleaned data are missing MEG channels: %s', strjoin(missing, ', '));
end
if numel(unique(cleanSelection)) ~= numel(cleanSelection)
    error('do_mne_ex_read_write_raw_SI_ver2:DuplicateChannelMap', ...
        'The cleaned MEG channel mapping contains duplicates.');
end

destinationChannels = megIndex(inputSelection);
sourceChannels = cleanSelection;
rawSampleCount = double(raw.last_samp - raw.first_samp + 1);
cleanSampleCount = size(clnData.trial{1}, 2);

if cleanSampleCount ~= rawSampleCount
    error('do_mne_ex_read_write_raw_SI_ver2:SampleCountMismatch', ...
        ['Cleaned and original recordings have different lengths: ' ...
         '%d versus %d samples.'], cleanSampleCount, rawSampleCount);
end
if abs(double(clnData.fsample)-double(raw.info.sfreq)) > ...
        10*eps(double(raw.info.sfreq))
    error('do_mne_ex_read_write_raw_SI_ver2:SamplingRateMismatch', ...
        ['Cleaned and original sampling rates differ: %.12g versus ' ...
         '%.12g Hz.'], clnData.fsample, raw.info.sfreq);
end

for k = 1:numel(megIndex)
    channel = megIndex(k);
    raw.info.chs(channel).cal = 1;
    raw.info.chs(channel).range = 1;
end

[outputFid, calibrations] = ...
    fiff_start_writing_raw(outfile, raw.info, picks);
outputOpen = true;
try
    calibrations(isMeg) = 1;
    blockSamples = max(1, round(double(cfg.blockSec)*raw.info.sfreq));
    nBlocks = ceil(rawSampleCount/blockSamples);
    fprintf('Streaming cleaned FIF in %d blocks of at most %.1f seconds.\n', ...
        nBlocks, blockSamples/raw.info.sfreq);

    firstBuffer = true;
    for block = 1:nBlocks
        cleanFirst = (block-1)*blockSamples + 1;
        cleanLast = min(block*blockSamples, rawSampleCount);
        rawFirst = double(raw.first_samp) + cleanFirst - 1;
        rawLast = double(raw.first_samp) + cleanLast - 1;

        [buffer, ~] = ...
            fiff_read_raw_segment(raw, rawFirst, rawLast, picks);
        buffer(destinationChannels, :) = ...
            clnData.trial{1}(sourceChannels, cleanFirst:cleanLast);

        if firstBuffer
            if rawFirst ~= 0
                fiff_write_int(outputFid, FIFF.FIFF_FIRST_SAMPLE, rawFirst);
            end
            firstBuffer = false;
        end

        fiff_write_raw_buffer(outputFid, buffer, calibrations);
        clear buffer
        if block == 1 || block == nBlocks || mod(block, 10) == 0
            fprintf('  export block %d/%d\n', block, nBlocks);
        end
    end

    fiff_finish_writing_raw(outputFid);
    outputOpen = false;
catch ME
    if outputOpen
        close_fif_safely(outputFid);
    end
    rethrow(ME);
end

clear inputCleanup
fprintf('Finished streaming cleaned FIF: %s\n', outfile);
end

function close_fif_safely(fid)
if isempty(fid) || ~isnumeric(fid) || fid < 0
    return;
end
try
    fclose(fid);
catch
end
end
