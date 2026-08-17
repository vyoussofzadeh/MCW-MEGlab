function do_mne_ex_read_write_raw_SI(cfg)
% do_mne_ex_read_write_raw_SI
% Write a FIF file using cleaned FieldTrip data (cfg.cln_data) while
% preventing mag/grad scale changes by forcing MEG calibration to unity.

infile   = cfg.infile;
outfile  = cfg.outfile;
cln_data = cfg.cln_data;

global FIFF;
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

% Read raw header/info
raw = fiff_setup_read_raw(infile);
picks = 1:length(raw.info.chs);

% Identify MEG channels
is_meg = false(size(picks));
for k = 1:numel(picks)
    is_meg(k) = (raw.info.chs(k).kind == FIFF.FIFFV_MEG_CH);
end
meg_idx = find(is_meg);

% Force MEG calibration to unity so output is interpreted as SI units
for ii = 1:numel(meg_idx)
    k = meg_idx(ii);
    raw.info.chs(k).cal   = 1;
    raw.info.chs(k).range = 1;
end

% Start writing
[outfid, cals] = fiff_start_writing_raw(outfile, raw.info, picks);
cals(is_meg) = 1;

% Read full segment and replace MEG samples with cleaned data
to = raw.last_samp;
first = 1;
last  = to;

[data, ~] = fiff_read_raw_segment(raw, first, last, picks);

% Map by labels (safer than assuming 1:306)
hdr = ft_read_header(infile);  % requires FieldTrip on path
[sel_in, sel_cln] = match_str(hdr.label(is_meg), cln_data.label);

data(meg_idx(sel_in), :) = cln_data.trial{1}(sel_cln, :);

% Write and finish
fiff_write_raw_buffer(outfid, data, cals);
fiff_finish_writing_raw(outfid);

% Close input file handle
% fclose(raw.fid);
end
