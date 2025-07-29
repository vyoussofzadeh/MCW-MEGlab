sStat   = stats_anim_pt;          % field that contains the structure

% --- Pick the field you want to average ---------------------------------
mapField = 'tmap';                      % 'tmap', 'pmap', or whichever you need

% --- Define the time window (in seconds, relative to file time zero) ----
t1 = 0.10;                              % start of window  (e.g., 100 ms)
t2 = 0.30;                              % end of window    (e.g., 300 ms)

t1 = 0.15;                              % start of window  (e.g., 100 ms)
t2 = 0.30;                              % end of window    (e.g., 300 ms)

t1 = 0.30;                              % start of window  (e.g., 100 ms)
t2 = 0.45;                              % end of window    (e.g., 300 ms)

t1 = 0.45;                              % start of window  (e.g., 100 ms)
t2 = 0.60;                              % end of window    (e.g., 300 ms)

t1 = 0.6;                              % start of window  (e.g., 100 ms)
t2 = 0.75;                              % end of window    (e.g., 300 ms)

t1 = 0.75;                              % start of window  (e.g., 100 ms)
t2 = 0.9;                              % end of window    (e.g., 300 ms)

t1 = 0.9;                              % start of window  (e.g., 100 ms)
t2 = 1.050;                              % end of window    (e.g., 300 ms)
        
t1 = 1.050; t2 = 1.2;

t1 = 1.2; t2 = 1.350;

t1 = 1.35; t2 = 1.5;

t1 = 1.5; t2 = 1.650;

idxTime = find(sStat.Time >= t1 & sStat.Time <= t2);   % logical index

% --- Average across those columns ---------------------------------------
avgMap      = mean(sStat.(mapField)(:,idxTime), 2);    % nVertices × 1
avgTime     = mean(sStat.Time(idxTime));               % single time stamp

%%

nTime   = numel(sStat.Time);            % 196 in your file
repTmap = repmat(avgMap, 1, nTime);     % nVertices × 196, every column identical

% ---- build a Brainstorm-compatible structure ---------------------------
sRep               = sStat;             % start from the original stats structure
sRep.tmap          = repTmap;           % replicated t-values
sRep.pmap          = sStat.pmap;% (optional) fill p-map the same size
sRep.df            = nan(size(repTmap));%           or replicate df likewise

% comment + bookkeeping
twin               = [t1,t2];       % <-- your averaging window (edit)
sRep.Comment       = sprintf('%s | mean %g %g s (replicated)', ...
                              sStat.Comment, twin(1), twin(2));
sRep.TimeBands     = twin;              % store the window you averaged