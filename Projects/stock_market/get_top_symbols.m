function [symbols, T] = get_top_symbols(indexName, saveCsv)
% GET_TOP_SYMBOLS  S&P 100 or Nasdaq-100 constituents from Wikipedia (no toolboxes).
% Usage:
%   [symbols, T] = get_top_symbols('sp100', 1);
%   [symbols, T] = get_top_symbols('nasdaq100', 1);

if nargin < 1 || isempty(indexName), indexName = 'sp100'; end
if nargin < 2, saveCsv = false; end
indexName = lower(string(indexName));

switch indexName
    case "sp100"
        url = 'https://en.wikipedia.org/wiki/S%26P_100';
        src = 'S&P 100 (Wikipedia)';
    case "nasdaq100"
        url = 'https://en.wikipedia.org/wiki/Nasdaq-100';
        src = 'Nasdaq-100 (Wikipedia)';
    otherwise
        error('indexName must be ''sp100'' or ''nasdaq100''.');
end

html = webread(url);                               % raw HTML
tblHtml = extract_wikitable_with_symbol(html);     % Primary approach
if isempty(tblHtml)
    % Try looser mode: split all wikitables and pass each to readtable via a temp file
    candidates = regexp(html, '(?is)<table[^>]*class="[^"]*wikitable[^"]*"[^>]*>.*?</table>', 'match');
    [T, ok] = try_each_table_with_readtable(candidates);
    if ~ok
        if indexName == "sp100"
            % Fallback: try Nasdaq-100 automatically (often parses reliably)
            warning('Could not parse S&P 100 table. Falling back to Nasdaq-100.');
            [symbols, T] = get_top_symbols('nasdaq100', saveCsv);
            return
        else
            error('Could not locate a usable constituents table on: %s', url);
        end
    end
else
    % Parse the HTML of the chosen table (no Text Analytics Toolbox)
    T = parse_symbol_name_table(tblHtml);
end

% Normalize and cap to 100 rows
if ~ismember('Symbol', T.Properties.VariableNames)
    error('Parsed table lacks a "Symbol" column.');
end
T.Symbol = strtrim(T.Symbol);
T = T(~cellfun(@isempty,T.Symbol),:);
if height(T) > 100, T = T(1:100,:); end

% Ensure Company column exists
if ~ismember('Company', T.Properties.VariableNames)
    if ismember('Security', T.Properties.VariableNames)
        T.Company = T.Security;
        T = movevars(T, 'Company', 'After', 'Symbol');
        T.Security = [];
    elseif ismember('Name', T.Properties.VariableNames)
        T.Company = T.Name;
        T = movevars(T, 'Company', 'After', 'Symbol');
        T.Name = [];
    else
        T.Company = repmat({''}, height(T), 1);
        T = movevars(T, 'Company', 'After', 'Symbol');
    end
end

% Add metadata
T.Source      = repmat({char(src)}, height(T), 1);
T.RetrievedAt = repmat({datestr(now,'yyyy-mm-dd HH:MM:SS')}, height(T), 1);

symbols = T.Symbol;

if saveCsv
    outName = sprintf('top_symbols_%s_%s.csv', char(indexName), datestr(now,'yyyymmdd'));
    writetable(T, outName);
    fprintf('Saved: %s\n', outName);
end
end

% -------- helpers --------
function tblHtml = extract_wikitable_with_symbol(html)
% Pick the first wikitable whose header contains 'symbol' or 'ticker'
tables = regexp(html, '(?is)<table[^>]*class="[^"]*wikitable[^"]*"[^>]*>.*?</table>', 'match');
tblHtml = '';
for i = 1:numel(tables)
    hdrRow = regexp(tables{i}, '(?is)<tr\b[^>]*>.*?</tr>', 'match');
    if isempty(hdrRow), continue; end
    ths = regexp(hdrRow{1}, '(?is)<th\b[^>]*>(.*?)</th>', 'tokens');
    if isempty(ths), continue; end
    hdrCells = strip_tags(join_tokens(ths));
    hdrCellsNorm = lower(regexprep(hdrCells, '\s+', ' '));
    if any(contains(hdrCellsNorm, 'symbol')) || any(contains(hdrCellsNorm, 'ticker'))
        tblHtml = tables{i};
        return
    end
end
end

function T = parse_symbol_name_table(tblHtml)
% Parse a single HTML table, extracting Symbol and Company/Security/Name
rows = regexp(tblHtml, '(?is)<tr\b[^>]*>.*?</tr>', 'match');
if numel(rows) < 2, error('No data rows.'); end

% Header
ths = regexp(rows{1}, '(?is)<th\b[^>]*>(.*?)</th>', 'tokens');
hdr = lower(regexprep(strip_tags(join_tokens(ths)), '\s+',' '));
symCol  = find(contains(hdr,'symbol')|contains(hdr,'ticker'), 1, 'first');
nameCol = find(contains(hdr,'company')|contains(hdr,'security')|contains(hdr,'name'), 1, 'first');
if isempty(symCol), symCol = 1; end
if isempty(nameCol), nameCol = symCol + 1; end

S = {}; C = {};
for r = 2:numel(rows)
    tds = regexp(rows{r}, '(?is)<td\b[^>]*>(.*?)</td>', 'tokens');
    if isempty(tds), continue; end
    cells = strip_tags(join_tokens(tds));
    s = ""; if symCol <= numel(cells), s = cells{symCol}; end
    s = clean_symbol(s);
    if strlength(s) == 0
        a = regexp(rows{r}, '(?is)<a\b[^>]*>(.*?)</a>', 'tokens');
        if ~isempty(a), s = clean_symbol(strip_tags(join_tokens({a{1}}))); end
    end
    if strlength(s) == 0, continue; end
    S{end+1,1} = char(s); %#ok<AGROW>

    nm = "";
    if nameCol <= numel(cells), nm = strtrim(cells{nameCol}); end
    C{end+1,1} = char(nm); %#ok<AGROW>
end

T = table(S, C, 'VariableNames', {'Symbol','Company'});
end

function [T, ok] = try_each_table_with_readtable(tableBlocks)
% Try to feed each <table>...</table> to readtable via a temp file.
ok = false; T = table();
for i = 1:numel(tableBlocks)
    try
        tmp = [tempname '.html'];
        fid = fopen(tmp,'w'); fwrite(fid, tableBlocks{i}); fclose(fid);
        Ti = readtable(tmp, 'TextType','string', 'PreserveVariableNames',true);
        delete(tmp);
        % Heuristic: keep table that has a 'Symbol' or 'Ticker' column
        cols = lower(string(Ti.Properties.VariableNames));
        if any(cols == "symbol") || any(contains(cols,"ticker"))
            % Normalize column names a little
            Ti.Properties.VariableNames = matlab.lang.makeValidName(Ti.Properties.VariableNames);
            % If variable names were row 1, fix by promoting first row
            if all(ismissing(Ti{1,:}))
                % do nothing
            end
            % Try to coerce Symbol/Company
            if ~ismember('Symbol', Ti.Properties.VariableNames)
                % look for something like 'Ticker_symbol' etc.
                idx = find(strcmpi(Ti.Properties.VariableNames, 'Ticker_symbol') | contains(cols, 'symbol') | contains(cols,'ticker'),1,'first');
                if ~isempty(idx), Ti.Symbol = string(Ti{:,idx}); end
            end
            if ~ismember('Company', Ti.Properties.VariableNames)
                cand = find(contains(cols, 'company')|contains(cols,'security')|contains(cols,'name'),1,'first');
                if ~isempty(cand), Ti.Company = string(Ti{:,cand}); end
            end
            if ismember('Symbol', Ti.Properties.VariableNames)
                T = table(cellstr(strtrim(string(Ti.Symbol))), cellstr(strtrim(string(Ti.Company))), ...
                          'VariableNames', {'Symbol','Company'});
                ok = true; return
            end
        end
    catch
        % try next
    end
end
end

function cells = join_tokens(tokCell)
cells = cellfun(@(t) strtrim(t{1}), tokCell, 'UniformOutput', false);
end

function s = strip_tags(x)
if iscell(x), s = cellfun(@strip_tags, x, 'UniformOutput', false); return; end
s = regexprep(x, '(?is)<script\b.*?</script>', '');
s = regexprep(s,  '(?is)<style\b.*?</style>',  '');
s = regexprep(s,  '(?is)<sup\b.*?</sup>',      '');
s = regexprep(s,  '<[^>]+>',                   '');
s = html_unescape(strtrim(s));
end

function s = html_unescape(s)
repl = {'&amp;','&'; '&nbsp;',' '; '&#160;',' '; '&ndash;','-'; '&mdash;','-'; '&middot;','·'};
for i = 1:size(repl,1), s = strrep(s, repl{i,1}, repl{i,2}); end
end

function s = clean_symbol(s)
s = regexprep(s, '\[\d+\]', '');
s = regexprep(s, '\s+',   '');
s = strrep(s, '','-'); s = strrep(s,'','-'); s = strrep(s, char(160),'');
end
