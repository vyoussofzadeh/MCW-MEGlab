function do_downsample_fif(varargin)
% batch_downsample_fif   down-sample all *Clean_raw.fif in tsss* subfolders
%
%  Usage examples
%     batch_downsample_fif                            % uses defaults
%     batch_downsample_fif('newFs', 200, 'force', true)
%
%  Namevalue pairs
%     'baseDir'   : root folder (default '/data/MEG/Research/SpikeDectection')
%     'pattern'   : wildcard for input files (default '*Clean_raw.fif')
%     'newFs'     : target sampling rate (Hz, default 200)
%     'outSuffix' : appended to basename (default '_DS')
%     'force'     : overwrite if outfile exists (default false)

% ------------------------------------------------- user params & defaults
p = inputParser;
p.addParameter('baseDir',   '/data/MEG/Research/SpikeDectection', @ischar);
p.addParameter('pattern',   '*Clean_raw.fif',                      @ischar);
p.addParameter('newFs',     200,                                   @isnumeric);
p.addParameter('outSuffix', '_DS',                                 @ischar);
p.addParameter('force',     false,                                 @islogical);
p.parse(varargin{:});
cfg = p.Results;

fprintf('[batch] Root = %s\n', cfg.baseDir);
fprintf('[batch] Pattern = %s   New Fs = %g Hz\n\n', ...
        cfg.pattern, cfg.newFs);

% ------------------------------------------------- discover input files
subDirs   = dir(fullfile(cfg.baseDir, 'tsss*'));
fileList  = {};

for d = 1:numel(subDirs)
    if ~subDirs(d).isdir, continue; end
    fifFiles = dir(fullfile(cfg.baseDir, subDirs(d).name, cfg.pattern));
    for f = 1:numel(fifFiles)
        fileList{end+1} = fullfile(fifFiles(f).folder, fifFiles(f).name); %#ok<AGROW>
    end
end

if isempty(fileList)
    fprintf('[batch] No files found.  Nothing to do.\n');
    return
end

fprintf('[batch] %d file(s) queued\n\n', numel(fileList));

% ------------------------------------------------- iterate
for idx = 1:numel(fileList)
    inFile  = fileList{idx};
    [folder, base, ext] = fileparts(inFile);                 % #ok<ASGLU>
    outFile = fullfile(folder, [base, cfg.outSuffix, ext]);

    if exist(outFile,'file') && ~cfg.force
        fprintf('[%02d/%02d] SKIP (exists)  %s\n', idx, numel(fileList), outFile);
        continue
    end

    tStart = tic;

    try
        % ---------- FieldTrip read + resample
        ft_defaults;         % ensure paths
        hdr = ft_read_header(inFile);   % quick header read
        oldFs = hdr.Fs;

        dcfg = []; dcfg.dataset = inFile; dcfg.demean = 'yes';
        data = ft_preprocessing(dcfg);

        rcfg = []; rcfg.resamplefs = cfg.newFs;
        data_ds = ft_resampledata(rcfg, data);

        % ---------- write out with FIFF I/O
        wcfg = [];                                  % for helper writer
        wcfg.infile    = inFile;
        wcfg.outfile   = outFile;
        wcfg.sfreq_new = cfg.newFs;
        wcfg.blocksec  = Inf;                       % process whole file
        wcfg.cln_data  = cell2mat(data_ds.trial);   % (nChan × nSamples)

        do_mne_write_downsample_whole(wcfg);

        % ---------- QC summary
        hdr_ds = ft_read_header(outFile);
        fprintf('[%02d/%02d] %s  |  %.0f ? %.0f Hz,  %d ? %d samp  (%.1f s)\n', ...
            idx, numel(fileList), base, ...
            oldFs, hdr_ds.Fs, hdr.nSamples, hdr_ds.nSamples, toc(tStart));

    catch ME
        warning('[%02d/%02d] FAILED %s  (%s)', idx, numel(fileList), inFile, ME.message);
        continue
    end
end

fprintf('\n[batch] Done.\n');
end
