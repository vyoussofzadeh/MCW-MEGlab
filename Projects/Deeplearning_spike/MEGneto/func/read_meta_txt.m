function meta = read_meta_txt(meta_fn, Fs, verbose)
%READ_META_TXT  Read reviewer meta_data.txt and return table + summary.
%
% Supports the common 3-column format:
%   sample_index   agree(0/1)   code(eg 5555)
%
% Inputs
%   meta_fn  : full path to *_meta_data.txt
%   Fs       : (optional) sampling rate for sample->seconds conversion
%   verbose  : (optional) true/false (default false)
%
% Output struct fields
%   meta.fn
%   meta.table    (table with columns sample, agree, code, [time_sec])
%   meta.header   (reserved; empty for this format)
%   meta.summary  (nRows, nApproved, nRejected, approvalRate)

    if nargin < 2 || isempty(Fs), Fs = NaN; end
    if nargin < 3 || isempty(verbose), verbose = false; end

    meta = struct();
    meta.fn = meta_fn;
    meta.header = struct();
    meta.table  = table();
    meta.summary = struct('nRows',0,'nApproved',NaN,'nRejected',NaN,'approvalRate',NaN);

    if ~isfile(meta_fn)
        if verbose
            fprintf('[read_meta_txt] Missing file: %s\n', meta_fn);
        end
        return;
    end

    % --- read file into lines (compatible with older MATLAB) ---
    txt = fileread(meta_fn);
    lines = string(regexp(txt, '\r\n|\n|\r', 'split'))';
    lines = lines(strlength(strtrim(lines)) > 0);

    if isempty(lines)
        if verbose
            fprintf('[read_meta_txt] Empty file: %s\n', meta_fn);
        end
        return;
    end

    % --- Try parse as strict 3-column numeric: sample agree code ---
    % Join all lines and sscanf into 3 x N, then transpose to N x 3
    M = sscanf(join(lines, newline), '%f %f %f', [3 Inf])';

    % Validate: each row should be integers-ish and agree should be 0/1
    ok = ~isempty(M) && size(M,2) == 3 && all(ismember(round(M(:,2)), [0 1]));

    if ok
        % Coerce to integers (these should be ints)
        sample = round(M(:,1));
        agree  = round(M(:,2));
        code   = round(M(:,3));

        Tmeta = table(sample, agree, code);

        % Add seconds if Fs provided
        if ~isnan(Fs) && Fs > 0
            Tmeta.time_sec = double(Tmeta.sample) ./ double(Fs);
        end

        meta.table = Tmeta;

        % Summary
        nCand     = height(Tmeta);
        nAgree    = sum(Tmeta.agree == 1);
        nDisagree = sum(Tmeta.agree == 0);

        meta.summary.nRows        = nCand;
        meta.summary.nApproved    = nAgree;
        meta.summary.nRejected    = nDisagree;
        meta.summary.approvalRate = nAgree / max(1, nCand);

        if verbose
            fprintf('[read_meta_txt] %s | N=%d agree=%d disagree=%d rate=%.3f\n', ...
                meta_fn, nCand, nAgree, nDisagree, meta.summary.approvalRate);
        end

        return;
    end

    % --- Fallback: store raw lines only (if format differs) ---
    meta.table = table(lines, 'VariableNames', {'raw'});
    meta.summary.nRows = height(meta.table);

    if verbose
        fprintf('[read_meta_txt] Unrecognized format, stored raw lines: %s\n', meta_fn);
    end
end
