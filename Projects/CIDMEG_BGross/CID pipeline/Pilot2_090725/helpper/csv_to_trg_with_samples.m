function csv_to_trg_with_samples(csvFile, outTrg, Fs, varargin)
p = inputParser;
addRequired(p,'csvFile',@(s)ischar(s)||isstring(s));
addRequired(p,'outTrg', @(s)ischar(s)||isstring(s));
addRequired(p,'Fs',     @(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'TimeUnit','s',@(s)any(strcmpi(s,{'s','ms'})));
addParameter(p,'T0',0,@(x)isnumeric(x)&&isscalar(x));
addParameter(p,'Tmax',[],@(x)isempty(x)||(isnumeric(x)&&isscalar(x)));
parse(p,csvFile,outTrg,Fs,varargin{:});
opt = p.Results;

T = readtable(csvFile,'VariableNamingRule','preserve');
labels = strtrim(string(T{:,1}));
times  = double(T{:,2});
if strcmpi(opt.TimeUnit,'ms'), times = times/1000; end

good   = labels~="" & ~isnan(times);
labels = labels(good); times = times(good);
[times, idx] = sort(times); labels = labels(idx);

samples = round((times - opt.T0) * opt.Fs);

if ~isempty(opt.Tmax)
    keep = (times >= opt.T0) & (times < opt.Tmax);
    dropped = sum(~keep);
    if dropped>0
        fprintf('Dropping %d events outside [%.3f, %.3f) s\n', dropped, opt.T0, opt.Tmax);
    end
    times = times(keep); labels = labels(keep); samples = samples(keep);
end

keep = samples >= 0;
if ~all(keep)
    fprintf('Dropping %d events with negative samples (time < T0)\n', sum(~keep));
    times = times(keep); labels = labels(keep); samples = samples(keep);
end

fid = fopen(char(outTrg),'w'); assert(fid>0,'Cannot write %s',outTrg);
fprintf(fid,'latency sample name\n');
for k=1:numel(times)
    fprintf(fid,'%.6f %d %s\n', times(k), samples(k), labels(k));
end
fclose(fid);
fprintf('Wrote %s with %d events.\n', outTrg, numel(times));
end
