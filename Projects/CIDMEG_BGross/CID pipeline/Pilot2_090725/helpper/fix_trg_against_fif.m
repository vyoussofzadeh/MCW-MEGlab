function fix_trg_against_fif(trgFile, fifFile, outTrg, varargin)
% fix_trg_against_fif  Align a .trg to a FIF: adjust for first_samp if needed and drop out-of-range events.
% Requires FieldTrip on path (ft_defaults).
%
% Usage:
%   fix_trg_against_fif('run1.trg','run1_raw.fif','run1_clamped.trg');
%   fix_trg_against_fif('run1.trg','run1_raw.fif','run1_clamped.trg','Verbose',true);

p = inputParser;
addRequired(p,'trgFile',@(s)ischar(s)||isstring(s));
addRequired(p,'fifFile',@(s)ischar(s)||isstring(s));
addRequired(p,'outTrg', @(s)ischar(s)||isstring(s));
addParameter(p,'Verbose',true,@islogical);
parse(p,trgFile,fifFile,outTrg,varargin{:});
opt = p.Results;

% --- Header (FieldTrip) ---
hdr = ft_read_header(char(fifFile));
Fs  = hdr.Fs;
nS  = hdr.nSamples;
firstSamp = 0;
if isfield(hdr,'FirstSample'), firstSamp = double(hdr.FirstSample); end
if isfield(hdr,'orig') && isfield(hdr.orig,'first_samp'), firstSamp = double(hdr.orig.first_samp); end
if isfield(hdr,'orig') && isfield(hdr.orig,'raw') && isfield(hdr.orig.raw,'first_samp')
    firstSamp = double(hdr.orig.raw.first_samp);
end
Tmax = (nS - 1) / Fs;   % run length in seconds, relative to file start (t=0)

% --- Read TRG ---
T = readtable(char(trgFile), FileType='text', Delimiter=' ', MultipleDelimsAsOne=true, ...
              ReadVariableNames=false, TextType='string');
if width(T) < 3
    error('TRG must have: latency sample name');
end
lat  = double(T.Var1);     % seconds (expected)
samp = double(T.Var2);     % sample indices
name = strtrim(string(T.Var3));

% --- Stats before ---
before_resp1 = sum(name=="resp1");
before_resp2 = sum(name=="resp2");

% --- Detect whether samp is absolute or relative ---
% If it's relative, samp  round(lat*Fs). If it's absolute, samp  firstSamp + round(lat*Fs).
rel_samp = round(lat * Fs);
abs_samp = firstSamp + rel_samp;

err_rel = median(abs(samp - rel_samp), 'omitnan');
err_abs = median(abs(samp - abs_samp), 'omitnan');

if err_abs < err_rel/2
    mode = "absolute"; new_samp = samp - firstSamp;   % convert to in-file samples
elseif err_rel < err_abs/2
    mode = "relative"; new_samp = samp;               % already in-file samples
else
    % ambiguous: default to relative if near-equal
    mode = "relative"; new_samp = samp;
end

% --- Convert to seconds from in-file samples to be consistent ---
lat2 = new_samp / Fs;

% --- Keep only events that land inside the file [0, nS-1] ---
keep = new_samp >= 0 & new_samp < nS;
dropped = find(~keep);

% --- Write clamped TRG (samples relative to file start; Brainstorm-friendly) ---
fid = fopen(char(outTrg),'w'); assert(fid>0, 'Cannot write %s', outTrg);
fprintf(fid,'latency sample name\n');
for k = 1:numel(lat2)
    if ~keep(k), continue; end
    fprintf(fid,'%.6f %d %s\n', lat2(k), new_samp(k), name(k));
end
fclose(fid);

% --- Report ---
after_resp1 = sum(name(keep)=="resp1");
after_resp2 = sum(name(keep)=="resp2");

if opt.Verbose
    fprintf('Fs=%.3f Hz, nSamples=%d, Tmax=%.3f s, firstSamp=%d\n', Fs, nS, Tmax, firstSamp);
    fprintf('TRG sample mode detected: %s (err_rel=%.1f, err_abs=%.1f)\n', mode, err_rel, err_abs);
    fprintf('Kept %d/%d events; dropped %d (outside file)\n', sum(keep), numel(keep), sum(~keep));
    fprintf('resp1: %d ? %d   resp2: %d ? %d\n', before_resp1, after_resp1, before_resp2, after_resp2);
    if ~isempty(dropped)
        fprintf('Dropped indices (1-based in TRG): '); fprintf('%d ', dropped); fprintf('\n');
    end
    fprintf('Wrote: %s\n', outTrg);
end
end