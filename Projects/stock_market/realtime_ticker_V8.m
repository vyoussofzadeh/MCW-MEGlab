function realtime_ticker_V8(symbols)
% Realtime Ticker (Twelve Data)  credit-safe, multi-symbol, round-robin
% Strict BUY/SELL:
%    uses last CLOSED bar (not the forming bar)
%    requires price confirmation (price > slow SMA for BUY; < for SELL)
% Extras:
%    5-minute trend confirmation (SMA10>SMA20 for BUY; reverse for SELL)
%    Cooldown after flips to reduce ping-pong
%    ATR(14) stop and 2R target on flips (status + table)
%    Scrollable table (fixed height), same width as chart (no horizontal scrollbar)

%% --- USER SETTINGS ---
% API_KEY  = '3debf9cf713846e8b240235b8536e294';
API_KEY  = '8bce5cd0881f426493f41753b3c244b6';
pollSec  = 8;                    % ~7.5 req/min (under 8/min cap)
batchN   = 1;                    % 1 symbol per tick (credit-safe)
interval = '1min';
fastN    = 20;
slowN    = 30;
cushion  = 0.001;                % 0.1%

% Strictness toggles (base)
USE_CLOSED_BAR           = true;   % decide on previous (closed) bar
REQUIRE_PRICE_ABOVE_SLOW = true;   % BUY needs price>slow SMA (SELL mirrors)

% Signal quality add-ons
USE_TF_CONFIRM   = true;        % require higher-TF trend confirm
TF_INTERVAL      = '5min';      % trend timeframe
TF_FAST          = 10;          % fast SMA on TF
TF_SLOW          = 20;          % slow SMA on TF

USE_COOLDOWN     = true;        % ignore re-flips for a bit after a flip
COOLDOWN_SEC     = 180;         % 3 minutes

USE_ATR_EXITS    = true;        % compute ATR stop/target on flips
ATR_LEN          = 14;
ATR_MULT         = 1.5;         % stop distance = ATR * 1.5
R_MULT           = 2.0;         % target = entry + 2R

% Table paging (visible rows per page)
PAGE_SIZE        = 50;

% ------------------------------------------------------------
if nargin < 1 || isempty(symbols)
    error('Pass a cell array of tickers, e.g., realtime_ticker_V8({''AAPL'',''MSFT''})');
end
symbols = cellfun(@(s) upper(strtrim(s)), symbols, 'uni', 0);

% Renderer (optional)
try, opengl hardware; catch, end
set(groot,'defaultFigureRenderer','painters');

% Kill any leftover timers
try, delete(timerfindall('Tag','realtime_ticker_timer')); end

% Quick key check (1 request)
if ~td_key_ok(API_KEY)
    error('Twelve Data API key invalid or blocked. Double-check it.');
end

%% --- UI (toolbar + table + chart) ---
f = uifigure('Name','Realtime Ticker (Twelve Data)','Position',[100 100 1200 680]);

% Root grid: 3 rows -> toolbar / table / chart
root = uigridlayout(f,[3 1]);
root.RowHeight   = {40, 340, '1x'};     % fixed toolbar + table height
root.ColumnWidth = {'1x'};
root.Padding       = [8 8 8 8];
root.RowSpacing    = 6;
root.ColumnSpacing = 6;

% Toolbar grid (dropdown | page text | status)
tb = uigridlayout(root,[1 3]);
tb.Layout.Row = 1; tb.Layout.Column = 1;
tb.ColumnWidth = {240, 160, '1x'};
tb.RowHeight   = {40};
tb.Padding       = [0 0 0 0];
tb.RowSpacing    = 0;
tb.ColumnSpacing = 8;

numAll   = numel(symbols);
numPages = max(1, ceil(numAll / PAGE_SIZE));
curPage  = 1;

[pageSyms, ~, ~] = getPageLists(symbols, curPage, PAGE_SIZE);
drop = uidropdown(tb,'Items',pageSyms,'Value',pageSyms{1});
drop.Layout.Row = 1; drop.Layout.Column = 1;

pageText = uilabel(tb,'Text',sprintf('Page %d / %d',curPage,numPages), ...
    'HorizontalAlignment','left');
pageText.Layout.Row = 1; pageText.Layout.Column = 2;

statusLbl = uilabel(tb,'Text','Initializing...');
statusLbl.Layout.Row = 1; statusLbl.Layout.Column = 3;

% Table spans full width of row 2 (vertical scroll enabled by fixed height)
colW = {45, 75, 65, 75, 65, 65, 70, 70, 70, 'auto'};
%        #   Sym Last Chg  S20  S50  Sig Stop Tgt  AsOf
uit = uitable(root, ...
    'ColumnName', {'#','Symbol','Last','Change %','SMA20*','SMA50*','Signal','Stop','Target','As of'}, ...
    'ColumnEditable', false(1,10), ...
    'RowStriping', 'on', ...
    'ColumnWidth', colW, ...
    'RowName', {});
uit.FontSize     = 11;
uit.Layout.Row   = 2;
uit.Layout.Column= 1;

% Chart spans full width of row 3
ax  = uiaxes(root); 
ax.Layout.Row    = 3; 
ax.Layout.Column = 1;
title(ax,'Price with SMA20/SMA50'); xlabel(ax,'Time'); ylabel(ax,'Price');

%% --- State ---
mk = @(s) matlab.lang.makeValidName(s);   % safe struct keys for tickers
symKeys = cellfun(mk, symbols, 'uni', 0);

dataCache    = struct;           % fields by key: ts, px, s1, s2, asOf, lastS1Val, lastS2Val
lastSignal   = containers.Map;   % key -> 'BUY'/'SELL'/'HOLD'
lastFlipTime = containers.Map;   % key -> datetime of last flip
lastPlan     = containers.Map;   % key -> struct('entry','stop','target','atr')
nextIdx      = 1;

% Credit window tracker
rate.windowStart = datetime('now'); rate.used = 0; rate.limit = 8; rate.safeMargin = 0;

% Initial fill
refreshFromCache();

%% --- Timer ---
t = timer('ExecutionMode','fixedSpacing','Period',pollSec, ...
          'TimerFcn',@tick,'ErrorFcn',@(~,e)disp(e),'Tag','realtime_ticker_timer');
cleanupObj = onCleanup(@cleanupEverything); %#ok<NASGU>
f.CloseRequestFcn = @onClose;
start(t);

%% ===================== NESTED HELPERS =====================

    function onClose(~,~)
        cleanupEverything();
        if isvalid(f), delete(f); end
    end

    function cleanupEverything()
        try, stop(t);   end
        try, delete(t); end
        try, delete(timerfindall('Tag','realtime_ticker_timer')); end
    end

    function refreshCreditWindow()
        if seconds(datetime('now') - rate.windowStart) >= 60
            rate.windowStart = datetime('now'); rate.used = 0;
        end
    end
    function tf = haveCredits(), tf = (rate.used + rate.safeMargin) < rate.limit; end

    function ensureCacheKey(key)
        if ~isfield(dataCache,key)
            dataCache.(key) = struct('ts',[],'px',[],'s1',[],'s2',[], ...
                                     'asOf',fmtAsOf(datetime('now')), ...
                                     'lastS1Val',NaN,'lastS2Val',NaN);
        end
    end

    function txt = asOfOrEmpty(key)
        if isfield(dataCache,key) && isfield(dataCache.(key),'asOf')
            txt = dataCache.(key).asOf;
        else
            txt = '';
        end
    end

    function [has, px] = cacheHasPx(key)
        has = isfield(dataCache,key) && isfield(dataCache.(key),'px') && ~isempty(dataCache.(key).px);
        px  = []; if has, px = dataCache.(key).px; end
    end

    function tf = cooldownActive(key)
        if ~USE_COOLDOWN || ~isKey(lastFlipTime, key), tf = false; return; end
        tf = seconds(datetime('now') - lastFlipTime(key)) < COOLDOWN_SEC;
    end

    function [upTrend, downTrend, tfAsOf, addCost] = tf_trend_confirm(sym)
        upTrend=false; downTrend=false; tfAsOf=fmtAsOf(datetime('now')); addCost=0;
        if ~USE_TF_CONFIRM, return; end
        [TTF, tfAsOf, ~] = td_getIntraday_small(sym, TF_INTERVAL, API_KEY);
        addCost = 1;
        if isempty(TTF) || height(TTF) < max(TF_FAST,TF_SLOW), return; end
        TTF = ensure_timestamp_var(TTF);
        pxTF = TTF.close;
        sF = movmean(pxTF, TF_FAST, 'Endpoints','shrink');
        sS = movmean(pxTF, TF_SLOW, 'Endpoints','shrink');
        upTrend   = sF(end) > sS(end);
        downTrend = sF(end) < sS(end);
    end

    function plan = atr_plan_from_bars(Tbars, idx, side)
        plan = struct('entry',NaN,'stop',NaN,'target',NaN,'atr',NaN);
        if ~USE_ATR_EXITS || height(Tbars) < max(ATR_LEN+1, idx), return; end
        Tbars = ensure_timestamp_var(Tbars);
        c = Tbars.close; h=Tbars.high; l=Tbars.low;
        prevC = [c(1); c(1:end-1)];
        TR = max([h-l, abs(h-prevC), abs(l-prevC)], [], 2);
        ATR = movmean(TR, ATR_LEN, 'Endpoints','shrink');
        a = ATR(idx); if ~(isfinite(a) && a>0), return; end
        entry = c(idx);
        if strcmp(side,'BUY')
            stop = entry - ATR_MULT*a;
            target = entry + R_MULT*(entry - stop);
        else
            stop = entry + ATR_MULT*a;
            target = entry - R_MULT*(stop - entry);
        end
        plan.entry=entry; plan.stop=stop; plan.target=target; plan.atr=a;
    end

    function [pageSymsOut, pageKeysOut, idxVec] = getPageLists(allSyms, page, pageSize)
        s = (page-1)*pageSize + 1;
        e = min(page*pageSize, numel(allSyms));
        idxVec = s:e;
        pageSymsOut = allSyms(idxVec);
        pageKeysOut = cellfun(@(x) matlab.lang.makeValidName(x), pageSymsOut, 'uni', 0);
    end

function refreshFromCache()
    [pageSymsLocal, pageKeysLocal] = getPageLists(symbols, curPage, PAGE_SIZE);
    drop.Items = pageSymsLocal; if ~isempty(pageSymsLocal), drop.Value = pageSymsLocal{1}; end
    pageText.Text = sprintf('Page %d / %d',curPage,numPages);

    % Table rows (10 columns now)
    Trows = cell(0,10);
    for k = 1:numel(pageSymsLocal)
        sym = pageSymsLocal{k}; key = pageKeysLocal{k};
        rowNum = (curPage-1)*PAGE_SIZE + k;  % absolute numbering across pages
        [hasPx, px] = cacheHasPx(key);

        if ~hasPx
            Trows(end+1,:) = {rowNum, sym,'n/a','n/a','n/a','n/a','n/a','n/a','n/a', asOfOrEmpty(key)}; %#ok<AGROW>
        else
            s1v = dataCache.(key).lastS1Val;
            s2v = dataCache.(key).lastS2Val;
            price  = px(end);
            change = (px(end)-px(1))/px(1)*100;
            [sigRow, ~, ~, ~, ~] = compute_signal_strict(px, fastN, slowN, cushion, USE_CLOSED_BAR, REQUIRE_PRICE_ABOVE_SLOW);

            stp = NaN; tgt = NaN;
            if isKey(lastPlan, key), p = lastPlan(key); stp = p.stop; tgt = p.target; end

            Trows(end+1,:) = {rowNum, sym, ...
                sprintf('%.2f',price), ...
                sprintf('%.2f',change), ...
                fmtOrNA(s1v), ...
                fmtOrNA(s2v), ...
                sigRow, ...
                fmtOrNA(stp), ...
                fmtOrNA(tgt), ...
                asOfOrEmpty(key)}; %#ok<AGROW>
        end
    end

        if isvalid(uit), uit.Data = Trows; end

        % Plot
        if isempty(pageSymsLocal), cla(ax); title(ax,'No symbols on this page'); return; end
        sel = drop.Value; cla(ax);
        iSel = find(strcmp(pageSymsLocal, sel), 1, 'first');
        if isempty(iSel), title(ax,'Select a symbol on this page'); return; end
        keySel = pageKeysLocal{iSel};
        [hasPxSel, pxSel] = cacheHasPx(keySel);
        if ~hasPxSel
            title(ax, sprintf('%s  no data yet', sel));
            text(ax, .5,.5,'No bars to plot yet','Units','normalized','HorizontalAlignment','center'); return
        end
        ts = dataCache.(keySel).ts; s1 = dataCache.(keySel).s1; s2 = dataCache.(keySel).s2;
        if numel(pxSel) < 2
            title(ax, sprintf('%s  waiting for more bars', sel));
            text(ax, .5,.5,'Waiting for more bars...','Units','normalized','HorizontalAlignment','center'); return
        end
        hold(ax,'on'); plot(ax, ts, pxSel, 'DisplayName','Price');
        if ~isempty(s1), plot(ax, ts, s1, 'DisplayName',sprintf('SMA%d',fastN)); end
        if ~isempty(s2), plot(ax, ts, s2, 'DisplayName',sprintf('SMA%d',slowN)); end
        hold(ax,'off'); legend(ax,'Location','best');
        title(ax, sprintf('%s  (%s)', sel, asOfOrEmpty(keySel)));
        if numel(ts)>1, ax.XLim = [ts(max(1,end-390)) ts(end)]; end

        statusLbl.Text = sprintf('Page %d/%d  credits used this minute: %d/%d', curPage, numPages, rate.used, rate.limit);
    end

    % ---- strict signal logic ----
    function [sig, s1, s2, lastS1Val, lastS2Val] = compute_signal_strict(px, nFast, nSlow, cush, useClosed, requirePriceConfirm)
        sig = 'HOLD'; lastS1Val = NaN; lastS2Val = NaN;
        if numel(px) < max(nFast,nSlow), s1=[]; s2=[]; return, end
        s1 = movmean(px, nFast, 'Endpoints','shrink');
        s2 = movmean(px, nSlow, 'Endpoints','shrink');

        k = numel(px); idx = k;
        if useClosed, if k<2, return; end, idx = k-1; end
        s1c = movmean(px(1:idx), nFast, 'Endpoints','shrink');
        s2c = movmean(px(1:idx), nSlow, 'Endpoints','shrink');
        lastS1Val = s1c(end); lastS2Val = s2c(end);
        pxDecide  = px(idx);

        buyCond  = ( lastS1Val > lastS2Val*(1+cush) );
        sellCond = ( lastS1Val < lastS2Val*(1-cush) );
        if requirePriceConfirm
            buyCond  = buyCond  && (pxDecide > lastS2Val);
            sellCond = sellCond && (pxDecide < lastS2Val);
        end
        if buyCond, sig='BUY'; elseif sellCond, sig='SELL'; else, sig='HOLD'; end
    end

    % ---- main timer tick ----
    function tick(~,~)
        try
            if ~isvalid(f), return; end
            refreshCreditWindow();

            % round-robin pick
            n = numel(symbols); if n==0, return; end
            m = min(batchN, n);
            fetchIdx = zeros(1,m);
            for j = 1:m
                fetchIdx(j) = nextIdx; nextIdx = nextIdx + 1; if nextIdx > n, nextIdx = 1; end
            end

            for ii = fetchIdx
                if ~haveCredits(), statusLbl.Text = 'Credit cap: pausing fetch until next minute...'; break; end
                sym  = symbols{ii};
                key  = symKeys{ii};

                [T, asOf, msg] = td_getIntraday_small(sym, interval, API_KEY);
                rate.used = rate.used + 1;

                if isempty(T) || height(T) < 2
                    ensureCacheKey(key); dataCache.(key).asOf = asOf;
                    if ~isempty(msg), statusLbl.Text = sprintf('%s: %s', sym, msg); end
                    continue
                end

                T = ensure_timestamp_var(T);

                % base signal
                px = T.close; ts = T.timestamp;
                [sig, s1, s2, s1v, s2v] = compute_signal_strict(px, fastN, slowN, cushion, USE_CLOSED_BAR, REQUIRE_PRICE_ABOVE_SLOW);

                % TF confirm (extra call)
                if USE_TF_CONFIRM && (strcmp(sig,'BUY') || strcmp(sig,'SELL'))
                    [upTF, downTF, ~, addCost] = tf_trend_confirm(sym);
                    rate.used = rate.used + addCost;
                    if strcmp(sig,'BUY')  && ~upTF,   sig='HOLD'; end
                    if strcmp(sig,'SELL') && ~downTF, sig='HOLD'; end
                end

                % cooldown
                if USE_COOLDOWN
                    isFlip = (~isKey(lastSignal,key) && ~strcmp(sig,'HOLD')) || ...
                             ( isKey(lastSignal,key) && ~strcmp(sig,lastSignal(key)) && ~strcmp(sig,'HOLD'));
                    if isFlip && cooldownActive(key)
                        if isKey(lastSignal,key), sig = lastSignal(key); else, sig = 'HOLD'; end
                    end
                end

                % cache
                dataCache.(key).ts   = ts;
                dataCache.(key).px   = px;
                dataCache.(key).s1   = s1;
                dataCache.(key).s2   = s2;
                dataCache.(key).asOf = asOf;
                dataCache.(key).lastS1Val = s1v;
                dataCache.(key).lastS2Val = s2v;

                % on flip -> plan + alert
                if ~isKey(lastSignal,key) || ~strcmp(lastSignal(key),sig)
                    lastSignal(key) = sig;
                    if ~strcmp(sig,'HOLD')
                        idx = numel(px); if USE_CLOSED_BAR, idx = max(1, idx-1); end
                        plan = atr_plan_from_bars(T, idx, sig);
                        lastFlipTime(key) = datetime('now');
                        lastPlan(key)     = plan;
                        try, beep; end %#ok<TRYNC>
                        if isfinite(plan.entry)
                            statusLbl.Text = sprintf('%s %s @ %.2f | stop %.2f | target %.2f | ATR %.3f', ...
                                sig, sym, plan.entry, plan.stop, plan.target, plan.atr);
                        else
                            statusLbl.Text = sprintf('%s %s', sig, sym);
                        end
                    end
                end
            end

            % Refresh current page only
            refreshFromCache();

        catch err
            disp('Timer tick error:'); disp(err.message);
        end
    end
end

%% ===================== SUBFUNCTIONS (shared) =====================

function [T, asOf, msg] = td_getIntraday_small(symbol, interval, API_KEY)
% Twelve Data fetch (~30 bars). Normalizes to 'timestamp' column.
msg = ''; T = table(); asOf = fmtAsOf(datetime('now'));
base = 'https://api.twelvedata.com/time_series';
opts = weboptions('Timeout',15);
try
    S = webread(base, 'symbol',symbol, 'interval',interval, ...
                      'outputsize','30', 'timezone','America/New_York', ...
                      'apikey',API_KEY, opts);
catch ME, msg = ['webread failed: ' ME.message]; return
end
if (isstruct(S) && isfield(S,'status') && strcmpi(S.status,'error')) || ...
   (isstruct(S) && isfield(S,'message') && ~isfield(S,'values'))
    if isfield(S,'message'), msg = char(S.message); else, msg = 'API error'; end, return
end
if ~isfield(S,'values') || isempty(S.values), msg='No values'; return, end

vals = S.values; n = numel(vals);
ts = NaT(n,1); op=zeros(n,1); hi=op; lo=op; cl=op; vo=op;
for i=1:n
    v=vals(i);
    ts(i)=datetime(v.datetime,'InputFormat','yyyy-MM-dd HH:mm:ss');
    op(i)=str2double(v.open); hi(i)=str2double(v.high); lo(i)=str2double(v.low);
    cl(i)=str2double(v.close); if isfield(v,'volume'), vo(i)=str2double(v.volume); else, vo(i)=NaN; end
end
T = table(ts, op, hi, lo, cl, vo, ...
          'VariableNames', {'timestamp','open','high','low','close','volume'});

T = ensure_timestamp_var(T);
T = sortrows(T,'timestamp');
asOf = fmtAsOf(T.timestamp(end));
end

function ok = td_key_ok(API_KEY)
ok = false; base = 'https://api.twelvedata.com/time_series';
opts = weboptions('Timeout',10);
try
    S = webread(base, 'symbol','AAPL', 'interval','1min', 'outputsize','1', ...
                      'timezone','America/New_York', 'apikey',API_KEY, opts);
    ok = isstruct(S) && isfield(S,'values') && ~isempty(S.values);
catch
    ok = false;
end
end

function T = ensure_timestamp_var(T)
% Normalize any plausible time column to 'timestamp'
if isempty(T), return; end
names = string(lower(T.Properties.VariableNames));
cand = find(ismember(names, ["timestamp","datetime","date","time","ts","var1"]), 1, 'first');
if isempty(cand)
    for i = 1:numel(names)
        try
            if isdatetime(T.(T.Properties.VariableNames{i}))
                cand = i; break
            end
        catch
        end
    end
end
if isempty(cand)
    error('Bars table has no recognizable time column. Columns: %s', strjoin(T.Properties.VariableNames, ', '));
end
if names(cand) ~= "timestamp"
    T.Properties.VariableNames{cand} = 'timestamp';
end
end

function s = fmtAsOf(ts)
% datestr tokens: month=mm (lower), minutes=MM (upper)
if ~isa(ts,'datetime'), ts = datetime(ts,'ConvertFrom','datenum'); end
s = datestr(ts, 'yyyy-mm-dd HH:MM');   % trimmed to fit
end

function txt = fmtOrNA(val)
if ~isfinite(val) || isnan(val)
    txt = 'n/a';
else
    txt = sprintf('%.2f', val);
end
end

