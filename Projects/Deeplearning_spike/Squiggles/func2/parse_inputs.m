function opts = parse_inputs(def, varargin)
% PARSE_INPUTS  Merge name/value pairs into an options struct.
% Usage:
%   opts = parse_inputs(defaultOpts, 'useGPU',true, 'norm','robust', ...)

p = inputParser;
p.CaseSensitive   = false;
p.PartialMatching = true;
p.KeepUnmatched   = true;   % keep unknowns (won't error; lets you add new fields later)

% Register all fields in the default struct so they're accepted as parameters
fn = fieldnames(def);
for i = 1:numel(fn)
    p.addParameter(fn{i}, def.(fn{i}));
end

% Parse incoming name/value pairs
p.parse(varargin{:});
opts = p.Results;

% Ensure boolean-like fields are logical
boolFields = intersect({'useGPU','savePlots','verbose','augment','useAug', ...
                        'usePCA','trainLSTM','trainCNN1D','trainRes1D','trainCNN2D','doFusion'}, fn);
for i = 1:numel(boolFields)
    b = opts.(boolFields{i});
    if isnumeric(b), b = logical(b); end
    opts.(boolFields{i}) = b;
end

% Basic sanitization
if isfield(opts,'bp') && ~isempty(opts.bp)
    bp = opts.bp(:).';
    assert(numel(bp)==2 && bp(1)>0, 'opts.bp must be [low high] in Hz');
    opts.bp = bp;
elseif isfield(opts,'bpHz') && ~isempty(opts.bpHz)
    bp = opts.bpHz(:).';
    assert(numel(bp)==2 && bp(1)>0, 'opts.bpHz must be [low high] in Hz');
    opts.bpHz = bp;
end

if isfield(opts,'notch') && (isempty(opts.notch) || opts.notch==0)
    opts.notch = [];
end
if isfield(opts,'notchHz') && (isempty(opts.notchHz) || opts.notchHz==0)
    opts.notchHz = [];
end

% Clamp percentiles
if isfield(opts,'tfixPct'), opts.tfixPct = max(50, min(100, round(opts.tfixPct))); end

% RNG seed (dont reset here; caller already set rng)
end
