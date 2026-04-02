function cleanup_MaxFilter_tSSS_forWrapper_ver2(subjIDpath, outDirpath)
% cleanup_MaxFilter_tSSS_forWrapper_ver2
% Headless wrapper: MaxFilter (t)SSS + SSP (ongoing/ECG/EOG) + write clean FIF
%
% subjIDpath: folder containing *raw.fif
% outDirpath: base output folder; products are written under outDirpath/sss/<runname>/

% Author: Vahab YoussofZadeh, 2026

% ---- normalize inputs
subjIDpath = char(subjIDpath);
outDirpath = char(outDirpath);

fprintf('Input file path: %s\n', subjIDpath);
fprintf('Output files will be saved in: %s\n', outDirpath);

if exist(subjIDpath,'dir') ~= 7
    error('Input folder not found: %s', subjIDpath);
end
if exist(outDirpath,'dir') ~= 7
    mkdir(outDirpath);
end

% --- Create meg_clinic XML templates (session.xml, raw.xml)
make_megclinic_xml_templates(subjIDpath);

% ---- analysis switches
do_analysis.ecg     = 1;
do_analysis.eog     = 0;
do_analysis.ongoing = 0;

% ---- Add paths (only if present)
pathList = { ...
    '/MEG_data/megclinic', ...
    '/MEG_data/MEG_Tools/Brainstorm/Brainstorm_2025/brainstorm3', ...
    '/MEG_data/MEG_Tools/MNE-2.7.0-3106-Linux-x86_64/share/matlab', ...
    '/MEG_data/MEG_Tools/mne', ...
    '/neuro/bin' ...
    };
for ii = 1:numel(pathList)
    if exist(pathList{ii}, 'dir') == 7
        addpath(genpath(pathList{ii}));
    else
        warning('Path not found: %s', pathList{ii});
    end
end

% ---- Work in input folder
cd(subjIDpath);
files = dir('*raw.fif');
if isempty(files)
    fprintf('No *raw.fif found in %s\n', subjIDpath);
    return;
end

for j = 1:numel(files)

    rawFile = files(j).name;
    baseRun = rawFile(1:end-8);   % remove '_raw.fif'
    baseStr = rawFile(1:end-4);   % remove '.fif'
    isEmpty = ~isempty(regexpi(baseStr, 'emptyroom'));

    % Output locations/names
    FN = make_fnames(outDirpath, baseRun, isEmpty);

    % Skip if final product exists (depends on empty vs non-empty)
    % Skip if final product exists (but still write XML if missing)
    analysesDir = fullfile(subjIDpath, 'analyses');
    rawNameInXml = rawFile;

    if isEmpty
        if exist(FN.sssFile,'file') == 2
            write_run_workflow_xml(FN.dir, baseRun, rawNameInXml, isEmpty, FN, analysesDir);
            fprintf('Skipping %s (SSS already exists).\n', baseStr);
            continue;
        end
    else
        if exist(fullfile(FN.dir, FN.cleanName),'file') == 2
            write_run_workflow_xml(FN.dir, baseRun, rawNameInXml, isEmpty, FN, analysesDir);
            fprintf('Skipping %s (already cleaned).\n', baseStr);
            continue;
        end
    end

    fprintf('\n=== Processing %s ===\n', baseStr);

    % 1) MaxFilter
    if isEmpty
        if exist(FN.sssFile, 'file') ~= 2
            runMaxFilterSSS(subjIDpath, rawFile, FN.sssFile);
        else
            fprintf('✓ SSS exists: %s\n', FN.sssFile);
        end

        % analysesDir = fullfile(subjIDpath, 'analyses');
        % rawNameInXml = rawFile;
        write_run_workflow_xml(FN.dir, baseRun, rawNameInXml, isEmpty, FN, analysesDir);

        fprintf('Empty-room: skipping SSP cleanup.\n');
        continue;
    else
        if exist(FN.tsssFile, 'file') ~= 2
            runMaxFilterTSSS(subjIDpath, rawFile, FN.tsssFile);
        else
            fprintf('✓ tSSS exists: %s\n', FN.tsssFile);
        end
    end

    % 2) Build proj files (must happen before applying)
    if do_analysis.ongoing
        ensure_ongoing_proj(FN.dir, FN.tsssName, FN.ongoingProjName);
    end
    if do_analysis.ecg
        ensure_ecg_proj(FN.dir, FN.tsssName, FN.ecgEveName, FN.ecgProjName);
    end
    if do_analysis.eog
        ensure_eog_proj(FN.dir, FN.tsssName, FN.eogEveName, FN.eogProjName);
    end

    % 3) Apply projections -> final ecgClean file (2 raw FIF outputs)
    appliedProj = {};
    if do_analysis.ongoing && exist(fullfile(FN.dir, FN.ongoingProjName), 'file') == 2
        appliedProj{end+1} = FN.ongoingProjName;
    end
    if do_analysis.ecg && exist(fullfile(FN.dir, FN.ecgProjName), 'file') == 2
        appliedProj{end+1} = FN.ecgProjName;
    end
    if do_analysis.eog && exist(fullfile(FN.dir, FN.eogProjName), 'file') == 2
        appliedProj{end+1} = FN.eogProjName;
    end

    write_clean_with_proj(FN.dir, FN.tsssName, FN.cleanName, appliedProj);

    % analysesDir = fullfile(subjIDpath, 'analyses');
    % rawNameInXml = rawFile;  % e.g. Run02_spont_supine_raw.fif
    write_run_workflow_xml(FN.dir, baseRun, rawNameInXml, isEmpty, FN, analysesDir);


    % fprintf('✓ Done: %s\n', fullfile(FN.dir, FN.cleanName));
    outPath = fullfile(FN.dir, FN.cleanName);
    if exist(outPath,'file')==2
        fprintf('✓ Done: %s\n', outPath);
    else
        fprintf('✗ Not created: %s\n', outPath);
    end

end

end

function FN = make_fnames(outRoot, baseRun, isEmpty)
FN = struct();

FN.dir = fullfile(outRoot, 'sss', baseRun);
if exist(FN.dir,'dir') ~= 7
    mkdir(FN.dir);
end

if isEmpty
    FN.sssName   = [baseRun '_raw_sss.fif'];
    FN.sssFile   = fullfile(FN.dir, FN.sssName);
    FN.cleanName = FN.sssName;   % <-- ADD THIS LINE
else
    FN.tsssName = [baseRun '_raw_t_sss.fif'];
    FN.tsssFile = fullfile(FN.dir, FN.tsssName);

    FN.cleanName       = [baseRun '_raw_t_sss_ecgClean_raw.fif'];
    FN.ecgEveName      = [baseRun '_raw_t_sss_ecg-eve.fif'];
    FN.ecgProjName     = [baseRun '_raw_t_sss_ecg-proj.fif'];

    FN.eogEveName      = [baseRun '_raw_t_sss_eog-eve.fif'];
    FN.eogProjName     = [baseRun '_raw_t_sss_eog-proj.fif'];

    FN.ongoingProjName = [baseRun '_raw_t_sss_ongoing-proj.fif'];
end
end

% ======================================================================
% MaxFilter runners (no -gui)
% ======================================================================
function runMaxFilterTSSS(inDir, rawFileName, outFileFull)
inputFile = fullfile(inDir, rawFileName);

fprintf('Running MaxFilter tSSS:\n  %s\n  -> %s\n', inputFile, outFileFull);

cmd = sprintf(['/neuro/bin/util/maxfilter -f %s -o %s ' ...
    '-ctc /neuro/databases/ctc/ct_sparse.fif ' ...
    '-cal /neuro/databases/sss/sss_cal.dat ' ...
    '-autobad off -st 10 -corr 0.9 -force'], inputFile, outFileFull);

[status, cmdout] = system(cmd);
if status ~= 0
    error('MaxFilter tSSS failed:\n%s', cmdout);
end
end

function runMaxFilterSSS(inDir, rawFileName, outFileFull)
inputFile = fullfile(inDir, rawFileName);

fprintf('Running MaxFilter SSS:\n  %s\n  -> %s\n', inputFile, outFileFull);

cmd = sprintf(['/neuro/bin/util/maxfilter -f %s -o %s ' ...
    '-ctc /neuro/databases/ctc/ct_sparse.fif ' ...
    '-cal /neuro/databases/sss/sss_cal.dat ' ...
    '-autobad off -force'], inputFile, outFileFull);

[status, cmdout] = system(cmd);
if status ~= 0
    error('MaxFilter SSS failed:\n%s', cmdout);
end
end

% ======================================================================
% SSP projection builders (MNE command line)
% ======================================================================
function ensure_ongoing_proj(runDir, rawName, projName)
projPath = fullfile(runDir, projName);
if exist(projPath,'file') == 2
    fprintf('✓ Ongoing proj exists: %s\n', projPath);
    return;
end

% Ongoing: whole-file projections
projTag = '_raw_sss_ongoing-proj';
cmd = sprintf(['mne_process_raw --cd %s --raw %s --makeproj --saveprojtag %s ' ...
    '--projnmag 3 --projngrad 3 --highpass 1.5 --lowpass 5 ' ...
    '--digtrigmask 0 --projgradrej -1 --projmagrej -1'], ...
    runDir, rawName, projTag);

run_unix_or_error(cmd, 'ONGOING makeproj failed');

% MNE writes *ongoing-proj.fif; rename newest one to our convention
rename_newest(runDir, '*_ongoing-proj.fif', projPath);
fprintf('✓ Ongoing proj created: %s\n', projPath);
end

function ensure_ecg_proj(runDir, rawName, eveName, projName)
evePath  = fullfile(runDir, eveName);
projPath = fullfile(runDir, projName);

% Event file
if exist(evePath,'file') ~= 2
    % We detect ECG events ourselves (no GUI). If no ECG channel -> skip.
    ok = write_ecg_events(runDir, rawName, evePath);
    if ~ok
        warning('ECG: No ECG events written; skipping ECG proj for %s', rawName);
        return;
    end
else
    fprintf('✓ ECG events exist: %s\n', evePath);
end

% Proj file
if exist(projPath,'file') == 2
    fprintf('✓ ECG proj exists: %s\n', projPath);
    return;
end

projTag = '_ecg-proj';
cmd = sprintf(['mne_process_raw --cd %s --raw %s --events %s --makeproj ' ...
    '--projtmin -0.08 --projtmax 0.08 --saveprojtag %s ' ...
    '--projnmag 1 --projngrad 1 --projevent 999 ' ...
    '--highpass 10 --lowpass 40 --digtrigmask 0'], ...
    runDir, rawName, eveName, projTag);

run_unix_or_error(cmd, 'ECG makeproj failed');
rename_newest(runDir, '*_ecg-proj.fif', projPath);
fprintf('✓ ECG proj created: %s\n', projPath);
end

function ensure_eog_proj(runDir, rawName, eveName, projName)
evePath  = fullfile(runDir, eveName);
projPath = fullfile(runDir, projName);

% If you have your own blinkDetect/get_eog pipeline you can drop it here.
% For now: if no eve exists, skip EOG safely.
if exist(evePath,'file') ~= 2
    warning('EOG: event file not present (%s). Skipping EOG SSP for %s', evePath, rawName);
    return;
end

if exist(projPath,'file') == 2
    fprintf('✓ EOG proj exists: %s\n', projPath);
    return;
end

% Determine event code in the EOG event file
events = mne_read_events(evePath);
eventNo = mode(double(events(:,3)));

projTag = '_eog-proj';
cmd = sprintf(['mne_process_raw --cd %s --raw %s --events %s --makeproj ' ...
    '--projtmin -0.2 --projtmax 0.2 --saveprojtag %s ' ...
    '--projnmag 1 --projngrad 1 --projevent %d ' ...
    '--highpass 1.5 --lowpass 15 --digtrigmask 0'], ...
    runDir, rawName, eveName, projTag, eventNo);

run_unix_or_error(cmd, 'EOG makeproj failed');
rename_newest(runDir, '*_eog-proj.fif', projPath);
fprintf('✓ EOG proj created: %s\n', projPath);
end

function write_clean_with_proj(runDir, rawName, outName, projList)

outPath = fullfile(runDir, outName);
if exist(outPath,'file') == 2
    fprintf('✓ Clean file exists: %s\n', outPath);
    return;
end
if isempty(projList)
    warning('No proj files to apply for %s (will not write clean file).', rawName);
    return;
end

stim = pick_exact_stim_name(fullfile(runDir, rawName));  % <-- NEW
mask = 0;

projStr = '';
for k = 1:numel(projList)
    projStr = [projStr ' --proj ' projList{k}];
end

cmd = sprintf(['mne_process_raw --cd %s --raw %s %s --projon --save %s ' ...
    '--digtrig "%s" --digtrigmask %d --filteroff'], ...
    runDir, rawName, projStr, outName, stim, mask);

run_unix_or_error(cmd, 'Applying proj + save clean failed');
end

function stim = pick_exact_stim_name(rawPath)
stim = 'STI101';  % preferred on your data
try
    f = fiff_setup_read_raw(rawPath);
    ch = string(f.info.ch_names);
    if any(ch == "STI101")
        stim = 'STI101';
    elseif any(ch == "STI 014")
        stim = 'STI 014';
    else
        % fallback: first STI* channel
        idx = find(startsWith(ch,"STI"), 1);
        if ~isempty(idx), stim = char(ch(idx)); end
    end
catch
end
end

% ======================================================================
% ECG event writing (batch-safe): detect ECG channel + QRS + write .fif events
% ======================================================================
function ok = write_ecg_events(runDir, rawName, evePath)
ok = false;

rawPath = fullfile(runDir, rawName);
try
    fiffsetup = fiff_setup_read_raw(rawPath);
catch ME
    warning('ECG: fiff_setup_read_raw failed: %s', ME.message);
    return;
end

[ecg, found] = get_ecg_channel_batchsafe(fiffsetup);
if ~found || isempty(ecg)
    warning('ECG: No ECG/EKG/BIO channel found in %s', rawName);
    return;
end

Options.sampRate = double(fiffsetup.info.sfreq);
Options.firstSamp = double(fiffsetup.first_samp);
Options.percentThresh = 60;
Options.noiseThresh = 2.5;
Options.maxCrossings = 3;
Options.minBeats = 10;
Options.ecgType = 999;

ecg_events = qrsDet_batchsafe(ecg, Options);
if isempty(ecg_events)
    warning('ECG: No QRS events detected in %s', rawName);
    return;
end

% Write MNE events (sample, 0, eventType)
firstSamp = Options.firstSamp;
adj = ecg_events(:) + firstSamp;

eventlist = zeros(numel(adj)+1,3);
eventlist(1,:) = [firstSamp 0 0];
eventlist(2:end,1) = adj;
eventlist(2:end,3) = Options.ecgType;

try
    mne_write_events(evePath, eventlist);
    ok = true;
catch ME
    warning('ECG: mne_write_events failed: %s', ME.message);
end
end

function [ecg, found] = get_ecg_channel_batchsafe(fiffsetup)
found = false;
ecg = [];

% ---- Robust channel selection (prefer ECG001 on new MEG)
chNames = strtrim(string(fiffsetup.info.ch_names));

ch = find(startsWith(upper(chNames), "ECG001"), 1);

if isempty(ch)
    cands = ["ECG","EKG","BIO001","BIO"];
    for i = 1:numel(cands)
        ch = find(contains(upper(chNames), upper(cands(i))), 1);
        if ~isempty(ch), break; end
    end
end

if isempty(ch)
    return;
end

start_samp = fiffsetup.first_samp;
end_samp   = fiffsetup.last_samp;
[ecg, ~] = fiff_read_raw_segment(fiffsetup, start_samp, end_samp, ch);

found = true;
end

function clean_events = qrsDet_batchsafe(ecg, Options)
clean_events = [];

minInterval   = round((60*Options.sampRate)/120);
BLANK_PERIOD  = minInterval;

ecg = ecg(:)';  % enforce 1xN

% bandpass via brainstorm (requires bst_bandpass_fft); if missing, fallback
if exist('bst_bandpass_fft','file') == 2
    filtecg = bst_bandpass_fft(ecg, Options.sampRate, 5, 35);
else
    % crude fallback: no filter
    filtecg = ecg;
end

absecg = abs(filtecg);
lenpts = length(absecg);

init = round(Options.sampRate);
if lenpts < 3*init
    return;
end

maxpt = [max(absecg(1:init)), max(absecg(init:2*init)), max(absecg(2*init:3*init))];
init_max = mean(maxpt);

thresh1 = init_max * (Options.percentThresh/100);
if thresh1 < 1e-6
    return;
end

k = 1; i = 1;
qrs_time = [];
qrs_rms = [];
qrs_numcross = [];

while i < lenpts - BLANK_PERIOD + 1
    if absecg(i) > thresh1
        window = absecg(i:i+BLANK_PERIOD);
        [~, maxTime] = max(window(1:floor(BLANK_PERIOD/2)));
        rmsv = sqrt(mean(window.^2));

        x = find(window > thresh1);
        y = diff(x);
        numcross = length(find(y > 1)) + 1;

        qrs_time(k) = maxTime + i; %#ok<AGROW>
        qrs_rms(k) = rmsv; %#ok<AGROW>
        qrs_numcross(k) = numcross; %#ok<AGROW>

        i = i + BLANK_PERIOD;
        k = k + 1;
    else
        i = i + 1;
    end
end

if isempty(qrs_time)
    return;
end

rms_mean   = mean(qrs_rms);
rms_std    = std(qrs_rms);
rms_thresh = rms_mean + (rms_std * Options.noiseThresh);

b = find(qrs_rms < rms_thresh);
a = qrs_numcross(b);
c = a < Options.maxCrossings;

clean_events = qrs_time(b(c));

if numel(clean_events) < Options.minBeats
    clean_events = [];
end
end

% ======================================================================
% Utilities
% ======================================================================
function run_unix_or_error(cmd, msg)
[status, out] = unix(cmd);
if status ~= 0
    error('%s:\n%s\nCMD:\n%s', msg, out, cmd);
end
end

function rename_newest(runDir, pattern, destFullPath)
d = dir(fullfile(runDir, pattern));
if isempty(d)
    return;
end
[~, idx] = max([d.datenum]);
src = fullfile(runDir, d(idx).name);
if ~strcmp(src, destFullPath)
    movefile(src, destFullPath);
end
end

function make_megclinic_xml_templates(runDir)
% Minimal XML templates in runDir (session.xml, raw.xml) like your previous code
analysesDir = fullfile(runDir, 'analyses');
reportDir   = fullfile(runDir, 'report');
if exist(analysesDir,'dir') ~= 7, mkdir(analysesDir); end
if exist(reportDir,'dir') ~= 7, mkdir(reportDir); end

sessionFile = fullfile(runDir, 'session.xml');
rawFile     = fullfile(runDir, 'raw.xml');

fid = fopen(sessionFile, 'w');
if fid < 0, error('Cannot write: %s', sessionFile); end
fprintf(fid, '<?xml version="1.0" encoding="utf-8"?>\n');
fprintf(fid, '<SESSION><TYPE>');
fprintf(fid, '<RAWDIR>%s</RAWDIR>', runDir);
fprintf(fid, '<MRIDIR/><EMPTY/><INTERICTALSPIKES/><INTERICTALFILES/><SPONT/>');
fprintf(fid, '<FUNCTIONALEVENTS/><FUNCTIONALFILES/><PROBEPOSITION/>');
fprintf(fid, '<BSTDIR/><BSTPROTOCOL/><BSTPROTOCOLTYPE/><BSTMRI/><BSTPIPELINE/>');
fprintf(fid, '<ANALYSIS>%s</ANALYSIS>', analysesDir);
fprintf(fid, '<REPORT>%s</REPORT>', reportDir);
fprintf(fid, '<FREQBANDS>delta: 1, 4; theta: 4, 8; alpha: 8, 12; beta: 13, 35; gamma1: 40, 80; gamma2: 80, 150; gamma3: 150, 250; gamma4: 250, 350; gamma5: 350, 550;</FREQBANDS>');
fprintf(fid, '<FXNLIMPORTSTRT/><FXNLIMPORTSTOP/><FXNLBASESTRT/><FXNLBASESTOP/>');
fprintf(fid, '</TYPE></SESSION>\n');
fclose(fid);

fid = fopen(rawFile, 'w');
if fid < 0, error('Cannot write: %s', rawFile); end
fprintf(fid, '<?xml version="1.0" encoding="utf-8"?>\n');
fprintf(fid, '<WORKFLOW><TYPE>');
fprintf(fid, '<RUNDIR>%s</RUNDIR>', runDir);
fprintf(fid, '<RAW>raw</RAW>');
fprintf(fid, '<AVE/><RAWSSS/><AVESSS/><ECGEVE/><ECGPROJ/><ECGCLEAN/>');
fprintf(fid, '<ARTEVE/><ARTPROJ/><ARTCLEAN/><EOGEVE/><EOGPROJ/><EOGCLEAN/>');
fprintf(fid, '<ONGNGPROJ/><ONGNGCLEAN/><EVENTLIST/><FXNLEVENTS/><AVEDESC/><AVECLEAN/>');
fprintf(fid, '<MASK/><STIMSOURCE/><USEREVE/><AVEUSEREVE/><CUSTOMAVEDESC/>');
fprintf(fid, '<MRIDIR/><MRIFIFF/><XFITDIR/><CFIT/><DIPOLES/><DIPOLEIMG/>');
fprintf(fid, '<BSTDIR/><BSTPROTOCOL/><BSTPROTOCOLTYPE/><BSTMRI/><BSTMEG/>');
fprintf(fid, '<REPORT/><XML/><EMTPY/><PROBEPOSITION/>');
fprintf(fid, '<ANALYSES>%s</ANALYSES>', analysesDir);
fprintf(fid, '</TYPE></WORKFLOW>\n');
fclose(fid);

fprintf('Wrote XML templates:\n  %s\n  %s\n', sessionFile, rawFile);
end


function write_run_workflow_xml(runDir, baseRun, rawName, isEmpty, FN, analysesDir)
% Writes <baseRun>.xml inside runDir, similar to your legacy per-run XML.

xmlName = [baseRun '.xml'];
xmlPath = fullfile(runDir, xmlName);

% Decide which files to reference
if isEmpty
    rawsss  = FN.sssName;             % e.g. Run02..._raw_sss.fif
    ecgeve  = [baseRun '_raw_sss_ecg-eve.fif'];
    ecgproj = [baseRun '_raw_sss_ecg-proj.fif'];
    ecgclean= [baseRun '_raw_sss_ecgClean_raw.fif'];
    fxnleve = [baseRun '_raw_sss-eve.fif'];   % optional, include if present
else
    rawsss  = FN.tsssName;            % for your wrapper this is tSSS file
    ecgeve  = FN.ecgEveName;
    ecgproj = FN.ecgProjName;
    ecgclean= FN.cleanName;
    fxnleve = '';                     % typically blank unless you generate one
end

% XFITDIR convention (as in your example)
xfitDir = fullfile(runDir, 'xfit');

fid = fopen(xmlPath, 'w');
if fid < 0
    warning('Could not write XML: %s', xmlPath);
    return;
end

fprintf(fid, '<?xml version="1.0" encoding="utf-8"?>');
fprintf(fid, '<WORKFLOW><TYPE>');

fprintf(fid, '<RUNDIR>%s</RUNDIR>', runDir);
fprintf(fid, '<RAW>%s</RAW>', rawName);

fprintf(fid, '<AVE/>');

% RAWSSS
if ~isempty(rawsss) && exist(fullfile(runDir, rawsss),'file')==2
    fprintf(fid, '<RAWSSS>%s</RAWSSS>', rawsss);
else
    fprintf(fid, '<RAWSSS/>');
end

fprintf(fid, '<AVESSS/>');

% ECG fields
if exist(fullfile(runDir, ecgeve),'file')==2
    fprintf(fid, '<ECGEVE>%s</ECGEVE>', ecgeve);
else
    fprintf(fid, '<ECGEVE/>');
end

if exist(fullfile(runDir, ecgproj),'file')==2
    fprintf(fid, '<ECGPROJ>%s</ECGPROJ>', ecgproj);
else
    fprintf(fid, '<ECGPROJ/>');
end

if exist(fullfile(runDir, ecgclean),'file')==2
    fprintf(fid, '<ECGCLEAN>%s</ECGCLEAN>', ecgclean);
else
    fprintf(fid, '<ECGCLEAN/>');
end

% Keep placeholders as in your example
fprintf(fid, '<ARTEVE/><ARTPROJ/><ARTCLEAN/>');
fprintf(fid, '<EOGEVE/><EOGPROJ/><EOGCLEAN/>');
fprintf(fid, '<ONGNGPROJ/><ONGNGCLEAN/>');
fprintf(fid, '<EVENTLIST/>');

% FXNLEVENTS (optional)
if ~isempty(fxnleve) && exist(fullfile(runDir, fxnleve),'file')==2
    fprintf(fid, '<FXNLEVENTS>%s</FXNLEVENTS>', fxnleve);
else
    fprintf(fid, '<FXNLEVENTS/>');
end

% AVEDESC (legacy)
fprintf(fid, '<AVEDESC>%s.ave</AVEDESC>', baseRun);
fprintf(fid, '<AVECLEAN/>');

% MASK/STIMSOURCE defaults
fprintf(fid, '<MASK>0</MASK>');
fprintf(fid, '<STIMSOURCE>0</STIMSOURCE>');

fprintf(fid, '<USEREVE/><AVEUSEREVE/><CUSTOMAVEDESC/>');
fprintf(fid, '<MRIDIR/><MRIFIFF/>');

fprintf(fid, '<XFITDIR>%s</XFITDIR>', xfitDir);
fprintf(fid, '<CFIT/><DIPOLES/><DIPOLEIMG/>');
fprintf(fid, '<BSTDIR/><BSTPROTOCOL/><BSTPROTOCOLTYPE/><BSTMRI/><BSTMEG/>');
fprintf(fid, '<REPORT/>');

% XML points to itself
fprintf(fid, '<XML>%s</XML>', xmlName);

fprintf(fid, '<EMTPY/><PROBEPOSITION/>');

% Analyses directory (you used dataset root /analyses)
fprintf(fid, '<ANALYSES>%s</ANALYSES>', analysesDir);

fprintf(fid, '</TYPE></WORKFLOW>');
fclose(fid);

fprintf('Wrote run XML: %s\n', xmlPath);
end
