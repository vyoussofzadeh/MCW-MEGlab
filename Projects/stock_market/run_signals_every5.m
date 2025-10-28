function run_signals_every5(symbols)
% RUN_SIGNALS_EVERY5  Console snapshot of signals every 5 minutes (no GUI).
% Usage:
%   run_signals_every5({'AAPL','MSFT','NVDA'})

% --- USER SETTINGS ---
API_KEY = '8bce5cd0881f426493f41753b3c244b6';
interval = '1min';
fastN    = 20;
slowN    = 30;
cushion  = 0.001;   % 0.1%

if nargin < 1 || isempty(symbols)
    error('Pass a cell array of tickers, e.g., run_signals_every5({''AAPL'',''MSFT''})');
end
symbols = upper(strtrim(symbols(:)'));  % row cell array

while true
    try
        T = snapshot_signals(symbols, interval, fastN, slowN, cushion, API_KEY);
        clc;
        fprintf('%s\n', datestr(now,'yyyy-mm-dd HH:MM:SS'));
        disp(T);
        % Optional logging:
        % writetable(T, 'signals_snapshot.csv', 'WriteMode','append');
    catch ME
        warning('Snapshot failed: %s', ME.message);
    end
    pause(5*60);   % 5 minutes
end
end

function T = snapshot_signals(symbols, interval, fastN, slowN, cushion, API_KEY)
% One-shot fetch + signal for a list of tickers (no timers, no GUI).
n = numel(symbols);
S = cell(n,1); Last = nan(n,1); ChgPct = nan(n,1);
SMAf = nan(n,1); SMAs = nan(n,1); Sig = cell(n,1); AsOf = cell(n,1);

for i = 1:n
    sym = symbols{i};
    [bars, asOf, msg] = td_getIntraday_small(sym, interval, API_KEY);
    if isempty(bars)
        S{i}   = sym;
        AsOf{i}= msg;
        Sig{i} = 'NODATA';
        continue
    end
    px = bars.close;
    S{i} = sym;
    AsOf{i} = asOf;
    Last(i) = px(end);
    ChgPct(i) = (px(end) - px(1)) / px(1) * 100;

    if numel(px) < max(fastN,slowN)
        Sig{i} = 'HOLD';
        continue
    end

    % Closed-bar decision + price confirmation
    k   = numel(px);
    idx = max(1, k-1);                 % last CLOSED bar
    s1c = movmean(px(1:idx), fastN, 'Endpoints','shrink');
    s2c = movmean(px(1:idx), slowN, 'Endpoints','shrink');
    smaF = s1c(end);  smaS = s2c(end);  p = px(idx);

    SMAf(i) = smaF;  SMAs(i) = smaS;

    sig = 'HOLD';
    if (smaF > smaS*(1+cushion)) && (p > smaS)
        sig = 'BUY';
    elseif (smaF < smaS*(1-cushion)) && (p < smaS)
        sig = 'SELL';
    end
    Sig{i} = sig;
end

T = table(S, Last, ChgPct, SMAf, SMAs, Sig, AsOf, ...
    'VariableNames', {'Symbol','Last','ChangePct','SMAfast','SMAslow','Signal','AsOf'});
end

function [T, asOf, msg] = td_getIntraday_small(symbol, interval, API_KEY)
% Twelve Data fetch (~30 bars). Toolbox-free, no string math.
msg  = '';
T    = table();
asOf = datestr(now,'yyyy-mm-dd HH:MM:SS');

base = 'https://api.twelvedata.com/time_series';
opts = weboptions('Timeout', 15);
try
    S = webread(base, 'symbol',symbol, 'interval',interval, ...
                      'outputsize','30', 'timezone','America/New_York', ...
                      'apikey',API_KEY, opts);
catch ME
    msg = ['webread failed: ' ME.message];
    return
end

% Handle API error / quota messages
if isstruct(S)
    hasValues = isfield(S,'values') && ~isempty(S.values);
    if ~hasValues
        if isfield(S,'message')
            msg = char(S.message);
        elseif isfield(S,'status')
            msg = ['API error: ' char(S.status)];
        else
            msg = 'API error / no values';
        end
        return
    end
else
    msg = 'Unexpected response';
    return
end

vals = S.values; n = numel(vals);
timestamp = NaT(n,1);
open = zeros(n,1); high = open; low = open; close = open; volume = open;

for i = 1:n
    v = vals(i);
    timestamp(i) = datetime(v.datetime,'InputFormat','yyyy-MM-dd HH:mm:ss');
    open(i)  = str2double(v.open);
    high(i)  = str2double(v.high);
    low(i)   = str2double(v.low);
    close(i) = str2double(v.close);
    if isfield(v,'volume'), volume(i) = str2double(v.volume); else, volume(i) = NaN; end
end

T = table(timestamp, open, high, low, close, volume, ...
          'VariableNames', {'timestamp','open','high','low','close','volume'});
T  = sortrows(T,'timestamp');
asOf = datestr(T.timestamp(end),'yyyy-mm-dd HH:MM:SS');
end
