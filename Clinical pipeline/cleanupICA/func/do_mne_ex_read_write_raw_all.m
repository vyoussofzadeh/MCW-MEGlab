function do_mne_ex_read_write_raw_all(cfg)
% do_mne_write_with_transform  —  copy a Neuromag FIF, adding dev→head-t
%
% cfg fields
%   .infile   : full path to original FIF (string)
%   .outfile  : desired output file  (string)
%   .T        : 4×4 rigid-body transform, *metres* (double)
%   .blocksec : (opt) seconds per read/write block, default 10
%
% Requires the MNE-Matlab FIFF I/O toolbox on the MATLAB path.
%
% V. Youssofzadeh & ChatGPT – 2025-06-30
% -------------------------------------------------------------------------

arguments
cfg struct
end
infile    = cfg.infile;
outfile   = cfg.outfile;
T         = cfg.T;                       % already metres – DON’T scale
blocksec  = getfield(cfg,'blocksec',10); % default 10 s

%% ------------------------------------------------------------------------
%  0)  Make sure FIFF constants are available
% -------------------------------------------------------------------------
global FIFF
if isempty(FIFF), FIFF = fiff_define_constants(); end
FDEV = int32(1);  FHEAD = int32(4);               % device / head codes

%% ------------------------------------------------------------------------
%  1)   Open the input file and patch the header
% -------------------------------------------------------------------------
raw = fiff_setup_read_raw(infile);

% ----- build dev_head_t struct ------------------------------------------
dev_head_t             = struct();
dev_head_t.from        = FDEV;
dev_head_t.to          = FHEAD;
dev_head_t.trans       = T;
dev_head_t.nelt        = int32(0);
dev_head_t.coord_id    = int32(0);

raw.info.dev_head_t    = dev_head_t;      % <<--- new transform

%% ------------------------------------------------------------------------
%  2)   Start writing the output file
% -------------------------------------------------------------------------
picks    = 1:numel(raw.info.chs);         % write all channels
[outfid,cals] = fiff_start_writing_raw(outfile, raw.info, picks);

first_sample    = raw.first_samp;
last_sample     = raw.last_samp;
block_samples   = blocksec * round(raw.info.sfreq);

fprintf('Copying %s  →  %s\n', infile, outfile);

% ----- stream through the file ------------------------------------------
for first = first_sample:block_samples:last_sample
    last = min(first + block_samples - 1, last_sample);
    
    [data, ~] = fiff_read_raw_segment(raw, first, last, picks);
    
    % -- write buffer -----------------------------------------------------
    if first == first_sample
        if first > 0
            fiff_write_int(outfid, FIFF.FIFF_FIRST_SAMPLE, first);
        end
    end
    fiff_write_raw_buffer(outfid, data, cals);
end

fiff_finish_writing_raw(outfid);
fclose(raw.fid);

fprintf('Finished.  dev_head_t embedded successfully.\n');
end


% function do_mne_ex_read_write_raw_all(cfg)
%
% infile = cfg.infile;
% outfile = cfg.outfile;
% cln_data = cfg.cln_data;
%
% %
% % function mne_ex_read_write_raw(infile,outfile);
% %
% % Read and write raw data in 60-sec blocks
% %
%
% %
% %   Author : Matti Hamalainen, MGH Martinos Center
% %   License : BSD 3-clause
% %
%
%
% global FIFF;
% if isempty(FIFF)
%    FIFF = fiff_define_constants();
% end
% %
% me = 'MNE:mne_ex_read_write_raw';
%
% try
%     raw = fiff_setup_read_raw(infile);
% catch
%     error(me,'%s',mne_omit_first_line(lasterr));
% end
%
%
% %
% picks = 1:length(raw.info.chs);
% [outfid,cals] = fiff_start_writing_raw(outfile,raw.info,picks);
% to          = raw.last_samp;
%
% %
% %   Read and write all the data
% first_buffer = true;
% for first = 1%from:quantum:to
%     last = to; %first+quantum-1;
%     if last > to
%         last = to;
%     end
%     try
%         [ data, ~ ] = fiff_read_raw_segment(raw,first,last,picks);
%     catch
%         fclose(raw.fid);
%         fclose(outfid);
%         error(me,'%s',mne_omit_first_line(lasterr));
%     end
%     %
%     %   You can add your own miracle here
%     %
%     data = cln_data.trial{1}; % replacing with raw data
%
%     fprintf(1,'Writing...');
%     if first_buffer
%        if first > 0
%            fiff_write_int(outfid,FIFF.FIFF_FIRST_SAMPLE,first);
%        end
%        first_buffer = false;
%     end
%     fiff_write_raw_buffer(outfid,data,cals);
%     fprintf(1,'[done]\n');
% end
%
% fiff_finish_writing_raw(outfid);
%
%
% %%%
%
%
% …
%
