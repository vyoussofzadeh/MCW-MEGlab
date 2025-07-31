function do_mne_write_downsample_clean(cfg)
%---------------------------------------------------------------------
% cfg fields
%   .infile      : source raw FIF
%   .outfile     : destination raw FIF
%   .sfreq_new   : target sampling rate (Hz)
%   .cln_data    : (opt) matrix to overwrite MEG channels (306 × nSamples)
%   .blocksec    : (opt) seconds per chunk  (default 15 s)
%---------------------------------------------------------------------

%% 0) Make sure FIFF constants are available  -------------------------
global FIFF
if isempty(FIFF)
    FIFF = fiff_define_constants();
end

if ~isfield(cfg,'blocksec'), cfg.blocksec = 15; end

%% 1) Open source file & patch header ---------------------------------
raw   = fiff_setup_read_raw(cfg.infile);
picks = 1:numel(raw.info.chs);

decim            = round(raw.info.sfreq / cfg.sfreq_new);
raw.info.sfreq   = cfg.sfreq_new;
raw.info.lowpass = cfg.sfreq_new/2;
raw.info.highpass= 0;

%% 2) Anti-alias FIR (design once)  -----------------------------------
fc   = 0.9*(cfg.sfreq_new/2);              % 90 % of new Nyquist
Wn   = fc / (raw.info.sfreq*decim/2);      % normalised
ord  = 10*decim;                           % taps  10×decim
b    = fir1(ord, Wn, hamming(ord+1));

%% 3) Start writer ----------------------------------------------------
[outfid,cals] = fiff_start_writing_raw(cfg.outfile, raw.info, picks);

blk_orig  = cfg.blocksec * decim * cfg.sfreq_new;   % samples in source
first_out = 0;
first_buffer = true;

for first = raw.first_samp : blk_orig : raw.last_samp
    last = min(first + blk_orig - 1, raw.last_samp);
    
    % -- read ---------------------------------------------------------
    data = fiff_read_raw_segment(raw, first, last, picks);
    
    % -- filter & decimate -------------------------------------------
    data_f  = filtfilt(b,1,data')';
    data_ds = data_f(:, 1:decim:end);
    
    % -- optional channel replacement --------------------------------
    if isfield(cfg,'cln_data') && ~isempty(cfg.cln_data)
        win_len = size(data_ds,2);
        data_ds(1:306,:) = cfg.cln_data(:, first_out+(1:win_len));
    end
    
    % -- write --------------------------------------------------------
    if first_buffer
        if first_out > 0
            fiff_write_int(outfid, FIFF.FIFF_FIRST_SAMPLE, first_out);
        end
        first_buffer = false;
    end
    fiff_write_raw_buffer(outfid, data_ds, cals);
    first_out = first_out + size(data_ds,2);
end

fiff_finish_writing_raw(outfid);
fclose(raw.fid);
fprintf('Done: %s ? %s  (%.0f Hz ? %.0f Hz)\n', ...
    cfg.infile, cfg.outfile, raw.info.sfreq*decim, cfg.sfreq_new);
end
