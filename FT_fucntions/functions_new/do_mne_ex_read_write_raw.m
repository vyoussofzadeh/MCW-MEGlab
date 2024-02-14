function do_mne_ex_read_write_raw(cfg)

infile = cfg.infile;
outfile = cfg.outfile;
cln_data = cfg.cln_data;

%
% function mne_ex_read_write_raw(infile,outfile);
%
% Read and write raw data in 60-sec blocks
%

%
%   Author : Matti Hamalainen, MGH Martinos Center
%   License : BSD 3-clause
%


global FIFF;
if isempty(FIFF)
   FIFF = fiff_define_constants();
end
%
me = 'MNE:mne_ex_read_write_raw';
%
% if nargin ~= 3
%     error(me,'Incorrect number of arguments');
% end
%
%   Setup for reading the raw data
%
try
    raw = fiff_setup_read_raw(infile);
catch
    error(me,'%s',mne_omit_first_line(lasterr));
end
%raw.info.projs = [];
%
%   Set up pick list: MEG + STI 014 - bad channels
%
%
% want_meg   = true;
% want_eog   = true;
% want_eeg   = true;
% want_stim  = true;
% include =  [ ];
% % include{1} = 'STI 014';
% 
% try
%     picks = fiff_pick_types(raw.info,want_meg,want_eeg,want_stim,raw.info.bads);
% catch
%     try
%         picks = fiff_pick_types(raw.info,want_meg,want_eeg,want_stim,include,raw.info.bads);
%     catch
%         error(me,'%s (channel list may need modification)',mne_omit_first_line(lasterr));
%     end
% end

%
picks = 1:length(raw.info.chs);
[outfid,cals] = fiff_start_writing_raw(outfile,raw.info,picks);
to          = raw.last_samp;

%
%   Read and write all the data
first_buffer = true;
for first = 1%from:quantum:to
    last = to; %first+quantum-1;
    if last > to
        last = to;
    end
    try
        [ data, ~ ] = fiff_read_raw_segment(raw,first,last,picks);
    catch
        fclose(raw.fid);
        fclose(outfid);
        error(me,'%s',mne_omit_first_line(lasterr));
    end
    %
    %   You can add your own miracle here
    %
    data1 = cln_data.trial{1}; % replacing with raw data
    data(1:306,:) = data1;
%     cals(1:306) = 1; 
    
    fprintf(1,'Writing...');
    if first_buffer
       if first > 0
           fiff_write_int(outfid,FIFF.FIFF_FIRST_SAMPLE,first);
       end
       first_buffer = false;
    end
    fiff_write_raw_buffer(outfid,data,cals);
    fprintf(1,'[done]\n');
end

fiff_finish_writing_raw(outfid);
% fclose(raw.fid);
