%% The CID (Pilot2) MEG

% Sanity check of trigger info - with manual inputs
% Writtern by MCW group, Youssofzadeh, Vahab <vyoussofzadeh@mcw.edu>
% Lastest update: 09/12/2025


%%
megdata = '/group/bgross/work/CIDMEG/ECOG_MEG_data/MEG/mcwa086_v1/250822/tsss/ICA/Run_02/ica_clean.fif';
megdata = '/group/bgross/work/CIDMEG/ECOG_MEG_data/MEG/mcwa086_v1/250822/tsss/ICA/Run_03/ica_clean.fif';
megdata = '/group/bgross/work/CIDMEG/ECOG_MEG_data/MEG/mcwa086_v1/250822/tsss/ICA/Run_04/ica_clean.fif';

% Header and events
hdr = ft_read_header(megdata);
evt = ft_read_event(megdata, 'detectflank','up');  % up/down/both as needed

% Keep only trigger-like events
types = {evt.type};
isTrig = ismember(types, {'trigger','STI101','STI001','STI 014','UPPT001','UPPT002'});
ev = evt(isTrig);

% Times in seconds
samples = double([ev.sample]);
firstSamp = 0;
if isfield(hdr,'FirstSample'), firstSamp = double(hdr.FirstSample);
elseif isfield(hdr,'orig') && isfield(hdr.orig,'first_samp'), firstSamp = double(hdr.orig.first_samp);
elseif isfield(hdr,'orig') && isfield(hdr.orig,'raw') && isfield(hdr.orig.raw,'first_samp')
    firstSamp = double(hdr.orig.raw.first_samp);
end
% t = (samples - firstSamp)/hdr.Fs;
t = (samples)/hdr.Fs;


vals = double([ev.value]);   % trigger codes

figure('Name','Trigger timeline'); 
stem(t, vals, 'filled');
xlabel('Time (s)'); ylabel('Trigger value'); grid on;
title('Triggers (value vs time)');

%%
evt = ft_read_event(megdata,'detectflank','up');
vals = double([evt.value]);
[u,~,ic] = unique(vals);
counts = accumarray(ic,1);
T = table(u',counts,'VariableNames',{'value','count'})


%%
% Requires FieldTrip on path (ft_defaults)
hdr   = ft_read_header(megdata);
Fs    = hdr.Fs;
Tmax  = (hdr.nSamples - 1) / Fs;    % end time in seconds
fprintf('Fs=%.3f Hz, nSamples=%d, Tmax=%.3f s\n', Fs, hdr.nSamples, Tmax);

%%
% Read your 2-col CSV: [label,time] (time in seconds or ms)
T = readtable('/group/bgross/work/CIDMEG/Analysis/Pilot2_090725/ResponseOnsets/cidmeg004_times_run1.csv', 'VariableNamingRule','preserve');
labels = string(T{:,1});
times  = double(T{:,2});   % if in ms: times = times/1000;

late = times >= Tmax;
fprintf('Events beyond raw end: %d\n', sum(late));
tab = table(labels(late), times(late), 'VariableNames', {'label','time_s'});
disp(tab)


%%
sum(T.Var1=="resp1"), sum(T.Var1=="resp2")
