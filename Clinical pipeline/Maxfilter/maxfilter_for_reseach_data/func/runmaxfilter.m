function Summary = runmaxfilter(sessionDir, baselineInput, bashScript, tsssDir, logDir)
% RUNMAXFILTER_VER2
% MaxFilter + QC + transform-report pipeline
%
% Features
% --------
% 1) Runs the bash MaxFilter wrapper
% 2) Saves bash output log
% 3) Creates headpos QC only when headpos txt files exist
% 4) Always creates transform reports from orig -> output FIF
% 5) Works with tsss / tsss_tans / tsss_tans_mvcomp
%
% Inputs
% ------
% sessionDir    : session folder containing orig/
% baselineInput : baseline FIF path or filename, or '' if not needed
% bashScript    : bash preprocessing wrapper
% tsssDir       : output folder (eg tsss, tsss_tans, tsss_tans_mvcomp)
% logDir        : folder for wrapper logs

    if nargin < 1 || isempty(sessionDir)
        sessionDir = pwd;
    end
    if nargin < 2
        baselineInput = '';
    end
    if nargin < 3 || isempty(bashScript)
        error('You must provide bashScript.');
    end
    if nargin < 4 || isempty(tsssDir)
        error('You must provide tsssDir.');
    end
    if nargin < 5 || isempty(logDir)
        logDir = fullfile(sessionDir, 'pipeline_logs');
    end

    if exist(sessionDir, 'dir') ~= 7
        error('Session directory not found: %s', sessionDir);
    end
    if exist(bashScript, 'file') ~= 2
        error('Bash script not found: %s', bashScript);
    end

    origDir = fullfile(sessionDir, 'orig');
    if exist(origDir, 'dir') ~= 7
        error('orig folder not found: %s', origDir);
    end

    if exist(logDir, 'dir') ~= 7
        mkdir(logDir);
    end

    % ---------------------------------------------------------------------
    % Resolve baseline
    % ---------------------------------------------------------------------
    baselineFile = '';
    baselineName = '';

    if ~isempty(baselineInput)
        baselineName = get_basename_only(baselineInput);

        if exist(baselineInput, 'file') == 2
            baselineFile = baselineInput;
        else
            candidate = fullfile(origDir, baselineName);
            if exist(candidate, 'file') == 2
                baselineFile = candidate;
            else
                baselineFile = '';
            end
        end
    end
    
    % ---------------------------------------------------------------------
    % Run bash wrapper
    % ---------------------------------------------------------------------
    [~, outFolderName] = fileparts(tsssDir);
    
    if ~isempty(baselineName)
        cmd = sprintf('cd %s && bash %s %s %s', ...
            shellquote(sessionDir), shellquote(bashScript), ...
            shellquote(baselineName), shellquote(outFolderName));
    else
        cmd = sprintf('cd %s && bash %s "" %s', ...
            shellquote(sessionDir), shellquote(bashScript), ...
            shellquote(outFolderName));
    end
    
%     % ---------------------------------------------------------------------
%     % Run bash wrapper
%     % ---------------------------------------------------------------------
%     if ~isempty(baselineName)
%         cmd = sprintf('cd %s && bash %s %s', ...
%             shellquote(sessionDir), shellquote(bashScript), shellquote(baselineName));
%     else
%         cmd = sprintf('cd %s && bash %s', ...
%             shellquote(sessionDir), shellquote(bashScript));
%     end
    
    fprintf('\nRunning MaxFilter bash wrapper...\n');
    fprintf('%s\n\n', cmd);
    
    
    [status, cmdout] = system(cmd);
    
    logFile = fullfile(logDir, 'Pipe_runmaxfilter.log');
    write_text_file(logFile, cmdout);
    
    fprintf('\nBash status: %d\n', status);
    fprintf('Bash output:\n%s\n', cmdout);
    
    if status ~= 0
        error('Bash wrapper failed with status %d.\nSee log: %s\n\nOutput:\n%s', status, logFile, cmdout);
    end
    
    if exist(tsssDir, 'dir') ~= 7
        error('Output folder not found after running MaxFilter: %s\nBash output was:\n%s', tsssDir, cmdout);
    end
    
%     [status, cmdout] = system(cmd);
%     
%     logFile = fullfile(logDir, 'Pipe_runmaxfilter.log');
%     write_text_file(logFile, cmdout);
%     
%     if status ~= 0
%         warning('Bash wrapper returned non-zero status: %d\nSee log: %s', status, logFile);
%     end
%     
%     if exist(tsssDir, 'dir') ~= 7
%         error('Output folder not found after running MaxFilter: %s', tsssDir);
%     end
    
    % ---------------------------------------------------------------------
    % QC output folder
    % ---------------------------------------------------------------------
    qcDir = fullfile(tsssDir, 'qc_jpeg');
    if exist(qcDir, 'dir') ~= 7
        mkdir(qcDir);
    end

    % ---------------------------------------------------------------------
    % Discover runs
    % ---------------------------------------------------------------------
    rawFiles = dir(fullfile(origDir, '*_raw.fif'));
    rawFiles = rawFiles(~[rawFiles.isdir]);

    if isempty(rawFiles)
        warning('No raw FIF files found in orig/: %s', origDir);
        Summary = [];
        return;
    end

    hpFiles = dir(fullfile(tsssDir, '*headpos*.txt'));
    hpFiles = hpFiles(~[hpFiles.isdir]);

    % map runStem -> headpos file
    hpMap = containers.Map();
    for k = 1:numel(hpFiles)
        hpPath = fullfile(hpFiles(k).folder, hpFiles(k).name);
        hpStem = local_stem(hpFiles(k).name);
        tok = regexp(hpStem, '^(.*)_headpos$', 'tokens');
        if ~isempty(tok)
            runStem = tok{1}{1};
        else
            runStem = hpStem;
        end
        hpMap(runStem) = hpPath;
    end

    if isempty(hpFiles)
        warning('No headpos txt files found in: %s', tsssDir);
    end

    % ---------------------------------------------------------------------
    % QC thresholds
    % ---------------------------------------------------------------------
    gofThr = 0.98;
    errThr = 0.005;   % meters = 5 mm

    % ---------------------------------------------------------------------
    % Per-run summary
    % ---------------------------------------------------------------------
    Summary = repmat(empty_summary_struct(), 0, 1);

    for k = 1:numel(rawFiles)
        runName = rawFiles(k).name;
        runStem = local_stem(runName);

        fprintf('[%d/%d] QC/Transform: %s\n', k, numel(rawFiles), runName);

        S = empty_summary_struct();
        S.run_name = runStem;

        % -------------------------------------------------------------
        % Headpos QC if available
        % -------------------------------------------------------------
        if isKey(hpMap, runStem)
            hpfile = hpMap(runStem);
            try
                S = qc_one_headpos_file(hpfile, qcDir, gofThr, errThr);
            catch ME
                warning('Headpos QC failed for %s\n%s', hpfile, ME.message);
                S = empty_summary_struct();
                S.hpfile = hpfile;
                S.run_name = runStem;
                S.note = append_note(S.note, ['Headpos QC failed: ' ME.message]);
            end
        else
            S.hpfile = '';
            S.note = append_note(S.note, 'No headpos txt found for this run');
        end

        % -------------------------------------------------------------
        % Transform report always
        % -------------------------------------------------------------
        try
            S = add_transform_report(S, sessionDir, tsssDir, baselineFile, qcDir);
        catch ME
            S.note = append_note(S.note, ['Transform export failed: ' ME.message]);
        end

        Summary(end+1,1) = S; %#ok<AGROW>
    end

    % ---------------------------------------------------------------------
    % Save summary
    % ---------------------------------------------------------------------
    save(fullfile(qcDir, 'QC_summary.mat'), 'Summary');
    write_summary_csv(fullfile(qcDir, 'QC_summary.csv'), Summary);

    fprintf('\nDone.\n');
    fprintf('QC JPEG folder: %s\n', qcDir);
    fprintf('Summary CSV   : %s\n', fullfile(qcDir, 'QC_summary.csv'));
    fprintf('Summary MAT   : %s\n', fullfile(qcDir, 'QC_summary.mat'));
end


function S = qc_one_headpos_file(hpfile, qcDir, gofThr, errThr)

    A = read_hp_txt_compat(hpfile);

    if isempty(A)
        error('Headpos file is empty or unreadable.');
    end

    A = A(~any(isnan(A),2), :);

    if isempty(A) || size(A,2) < 10
        error('Headpos file does not contain the expected 10 columns.');
    end

    t   = A(:,1);
    x   = A(:,5);
    y   = A(:,6);
    z   = A(:,7);
    gof = A(:,8);
    err = A(:,9);
    vel = A(:,10);

    if numel(t) > 1
        dt = diff(t);
        median_dt_s = median(dt);
    else
        median_dt_s = NaN;
    end

    bad = (gof < gofThr) | (err > errThr);
    d = sqrt((x - x(1)).^2 + (y - y(1)).^2 + (z - z(1)).^2);

    fname = local_stem(hpfile);
    jpgFile = fullfile(qcDir, [fname '_QC.jpg']);

    fig = figure('Visible', 'off', 'Color', 'w', 'Position', [100 100 900 700]);

    subplot(4,1,1);
    plot(t, [x y z], 'LineWidth', 1);
    ylabel('x / y / z');
    legend({'x','y','z'});
    title(fname, 'Interpreter', 'none');

    subplot(4,1,2);
    plot(t, gof, 'k', 'LineWidth', 1); hold on;
    plot([t(1) t(end)], [gofThr gofThr], 'r--', 'LineWidth', 1);
    ylabel('GOF');

    subplot(4,1,3);
    plot(t, err, 'k', 'LineWidth', 1); hold on;
    plot([t(1) t(end)], [errThr errThr], 'r--', 'LineWidth', 1);
    ylabel('Error');

    subplot(4,1,4);
    plot(t, vel, 'k', 'LineWidth', 1);
    ylabel('Velocity');
    xlabel('Time (s)');

    set(fig, 'PaperPositionMode', 'auto');
    print(fig, '-djpeg', '-r150', jpgFile);
    close(fig);

    S = empty_summary_struct();
    S.hpfile              = hpfile;
    S.run_name            = fname;
    S.rows                = size(A,1);
    S.median_dt_s         = median_dt_s;
    S.gof_median          = median(gof);
    S.gof_p01             = prctile(gof,1);
    S.error_median_mm     = 1000 * median(err);
    S.error_p99_mm        = 1000 * prctile(err,99);
    S.max_velocity        = max(vel);
    S.flagged_count       = sum(bad);
    S.flagged_percent     = 100 * mean(bad);
    S.max_displacement_mm = 1000 * max(d);
    S.jpg_file            = jpgFile;
    S.note                = '';
end


function S = empty_summary_struct()
    S = struct( ...
        'hpfile', '', ...
        'run_name', '', ...
        'rows', NaN, ...
        'median_dt_s', NaN, ...
        'gof_median', NaN, ...
        'gof_p01', NaN, ...
        'error_median_mm', NaN, ...
        'error_p99_mm', NaN, ...
        'max_velocity', NaN, ...
        'flagged_count', NaN, ...
        'flagged_percent', NaN, ...
        'max_displacement_mm', NaN, ...
        'jpg_file', '', ...
        'transform_txt', '', ...
        'note', '');
end


function write_summary_csv(csvFile, Summary)
    fid = fopen(csvFile, 'wt');
    if fid == -1
        error('Cannot write CSV file: %s', csvFile);
    end
    cleanupObj = onCleanup(@() fclose(fid)); 

    fprintf(fid, ['hpfile,run_name,rows,median_dt_s,gof_median,gof_p01,' ...
        'error_median_mm,error_p99_mm,max_velocity,flagged_count,' ...
        'flagged_percent,max_displacement_mm,jpg_file,transform_txt,note\n']);

    for i = 1:numel(Summary)
        s = Summary(i);

        fprintf(fid, '%s,%s,%.0f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.0f,%.3f,%.6f,%s,%s,%s\n', ...
            csv_escape(s.hpfile), ...
            csv_escape(s.run_name), ...
            s.rows, ...
            s.median_dt_s, ...
            s.gof_median, ...
            s.gof_p01, ...
            s.error_median_mm, ...
            s.error_p99_mm, ...
            s.max_velocity, ...
            s.flagged_count, ...
            s.flagged_percent, ...
            s.max_displacement_mm, ...
            csv_escape(s.jpg_file), ...
            csv_escape(s.transform_txt), ...
            csv_escape(s.note));
    end
end


function out = csv_escape(str)
    if isempty(str)
        out = '""';
        return;
    end
    str = strrep(str, '"', '""');
    out = ['"' str '"'];
end


function write_text_file(fname, txt)
    fid = fopen(fname, 'wt');
    if fid == -1
        warning('Could not write log file: %s', fname);
        return;
    end
    cleanupObj = onCleanup(@() fclose(fid)); 
    fprintf(fid, '%s', txt);
end


function nameOnly = get_basename_only(pathOrName)
    tok = regexp(pathOrName, '[^/\\]+$', 'match');
    if isempty(tok)
        nameOnly = pathOrName;
    else
        nameOnly = tok{1};
    end
end


function q = shellquote(s)
    s = strrep(s, '"', '\"');
    q = ['"' s '"'];
end


function d = local_dirname(p)
    if isempty(p)
        d = '';
        return;
    end

    idx = find((p == '/') | (p == '\'), 1, 'last');
    if isempty(idx)
        d = '';
    else
        d = p(1:idx-1);
    end
end


function s = local_stem(p)
    if isempty(p)
        s = '';
        return;
    end

    idx = find((p == '/') | (p == '\'), 1, 'last');
    if isempty(idx)
        name = p;
    else
        name = p(idx+1:end);
    end

    dotIdx = find(name == '.', 1, 'last');
    if isempty(dotIdx)
        s = name;
    else
        s = name(1:dotIdx-1);
    end
end


function A = read_hp_txt_compat(hpfile)
% Compatibility reader for MaxFilter head-position TXT files

    if exist('readmatrix', 'file') == 2
        A = readmatrix(hpfile, 'FileType', 'text');
        A = A(~any(isnan(A), 2), :);
        return;
    end

    fid = fopen(hpfile, 'rt');
    if fid == -1
        error('Could not open file: %s', hpfile);
    end
    cleaner = onCleanup(@() fclose(fid)); 

    A = [];
    while ~feof(fid)
        line = fgetl(fid);
        if ~ischar(line)
            continue;
        end

        line = strtrim(line);
        if isempty(line)
            continue;
        end

        vals = sscanf(line, '%f');

        if numel(vals) >= 10
            A(end+1, :) = vals(1:10).'; %#ok<AGROW>
        end
    end

    if isempty(A)
        warning('No numeric rows were read from: %s', hpfile);
    end
end


function S = add_transform_report(S, sessionDir, tsssDir, baselineFile, qcDir)

    runStem = S.run_name;
    tok = regexp(runStem, '^(.*)_headpos$', 'tokens');
    if ~isempty(tok)
        runStem = tok{1}{1};
    end

    origDir = fullfile(sessionDir, 'orig');

    beforeFile = find_matching_fif(origDir, runStem);
    afterFile  = find_matching_fif(tsssDir, runStem);

    S.transform_txt = '';

    if isempty(beforeFile)
        S.note = append_note(S.note, 'Transform export skipped: missing original FIF');
        return;
    end

    if isempty(afterFile)
        S.note = append_note(S.note, 'Transform export skipped: missing output FIF');
        return;
    end

    reportTxt = fullfile(qcDir, [runStem '_transform.txt']);

    if isempty(baselineFile) || exist(baselineFile, 'file') ~= 2
        baselineFile = '';
    end

    [ok, msg] = write_showfiff_transform_report(beforeFile, afterFile, baselineFile, reportTxt);

    if ok
        S.transform_txt = reportTxt;
        if ~isempty(msg)
            S.note = append_note(S.note, msg);
        end
    else
        S.note = append_note(S.note, ['Transform export failed: ' msg]);
    end
end


function fifFile = find_matching_fif(folder, stem)

    fifFile = '';

    if exist(folder, 'dir') ~= 7
        return;
    end

    d = dir(fullfile(folder, [stem '.fif']));
    d = d(~[d.isdir]);
    if ~isempty(d)
        fifFile = fullfile(folder, d(1).name);
        return;
    end

    d = dir(fullfile(folder, [stem '.fif.gz']));
    d = d(~[d.isdir]);
    if ~isempty(d)
        fifFile = fullfile(folder, d(1).name);
        return;
    end

    d = dir(fullfile(folder, [stem '*.fif*']));
    d = d(~[d.isdir]);
    if ~isempty(d)
        fifFile = fullfile(folder, d(1).name);
    end
end


function [ok, msg] = write_showfiff_transform_report(beforeFile, afterFile, baselineFile, outTxt)

    ok = false;
    msg = '';

    showFiffExe = find_showfiff();
    if isempty(showFiffExe)
        msg = 'show_fiff not found';
        return;
    end

    cmd1 = sprintf('%s -vt 222 %s', shellquote(showFiffExe), shellquote(beforeFile));
    cmd2 = sprintf('%s -vt 222 %s', shellquote(showFiffExe), shellquote(afterFile));

    [st1, txt1] = system(cmd1);
    [st2, txt2] = system(cmd2);

    txt3 = '';
    st3 = 0;
    if ~isempty(baselineFile)
        cmd3 = sprintf('%s -vt 222 %s', shellquote(showFiffExe), shellquote(baselineFile));
        [st3, txt3] = system(cmd3);
    end

    fid = fopen(outTxt, 'wt');
    if fid == -1
        msg = ['Cannot open output file: ' outTxt];
        return;
    end
    cleanupObj = onCleanup(@() fclose(fid)); 

    fprintf(fid, '=== ORIGINAL FILE (orig) ===\n%s\n\n', beforeFile);
    fprintf(fid, '%s\n', txt1);

    fprintf(fid, '\n\n=== PROCESSED FILE ===\n%s\n\n', afterFile);
    fprintf(fid, '%s\n', txt2);

    if ~isempty(baselineFile)
        fprintf(fid, '\n\n=== TRANSFORM REFERENCE FILE ===\n%s\n\n', baselineFile);
        fprintf(fid, '%s\n', txt3);
    end

    ok = true;

    if st1 ~= 0 || st2 ~= 0 || st3 ~= 0
        msg = 'show_fiff returned non-zero status for at least one file';
    end
end


function exe = find_showfiff()

    exe = '';

    candidates = { ...
        '/opt/neuromag/bin/util/show_fiff', ...
        '/neuro/bin/util/show_fiff'};

    for i = 1:numel(candidates)
        if exist(candidates{i}, 'file') == 2
            exe = candidates{i};
            return;
        end
    end
end


function out = append_note(oldNote, newNote)

    if isempty(newNote)
        out = oldNote;
        return;
    end

    if isempty(oldNote)
        out = newNote;
    else
        out = [oldNote ' | ' newNote];
    end
end