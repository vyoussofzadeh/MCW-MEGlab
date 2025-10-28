clc
clear
clear
clear functions
clear realtime_ticker

%%
[syms_sp, T] = get_top_symbols('sp100', 0);
fid = fopen('sp100_symbols.txt','w'); fprintf(fid, '%s\n', syms_sp{:}); fclose(fid);

[syms_nas, T] = get_top_symbols('nasdaq100', 1);   % writes CSV
fid = fopen('nasdaq100_symbols.txt','w'); fprintf(fid, '%s\n', syms_nas{:}); fclose(fid);

%% Ver old
% realtime_ticker_V6(syms_sp(1:25))
% realtime_ticker_V6(syms_sp(25:50))
% realtime_ticker_V6(syms_sp(51:100))
% 
% realtime_ticker_V6(syms_nas(1:10))

%% Ver new (tales all the credits)
% realtime_ticker_V8(syms_sp(1:11))
realtime_ticker_V8(syms_sp(1:50))
realtime_ticker_V8(syms_sp(51:100))

%%
owned_symbols = {'TSLA', ...      
           'CB',...
           'INTU',...
           'ABBV',...
           'GILD',...
           'META',...
           'AMGN'};  % Sector ETFs (Tech, Financials, Energy)

realtime_ticker_V8(owned_symbols)

