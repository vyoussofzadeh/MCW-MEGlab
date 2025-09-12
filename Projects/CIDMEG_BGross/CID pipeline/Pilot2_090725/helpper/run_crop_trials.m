% Assumes 'neuromag_meg' is already in your workspace (with .Channel)
ch = neuromag_meg.Channel;

% Identify MEG channels (keep MAG/GRAD, drop EEG/EOG/ECG/EMG/MISC/STIM/REF/COMB)
types = string({ch.Type});
isMEG = startsWith(types,"MEG") ...
        & ~startsWith(types,"MEG REF") ...
        & ~startsWith(types,"MEG COMB") ...
        & ~contains(types,["EEG","EOG","ECG","EMG","MISC","STIM"]);

iMEG = find(isMEG);
fprintf('Keeping %d MEG channels (expected ~306)\n', numel(iMEG));

% Edit the in-memory struct
neuromag_meg.Channel = ch(iMEG);
neuromag_meg.Comment = sprintf('Neuromag channels (%d) - MEG only', numel(iMEG));

% Save a Brainstorm-compatible channel file
% (Brainstorm is happy with a file that just contains variable 'Channel')
Channel = neuromag_meg.Channel; % #ok<NASGU>
% save('channel_vectorview306_acc1.mat','Channel','-v7');

% Optional: keep a full-style file (same fields as your original) instead:
% save('channel_vectorview306_onlyMEG_full.mat','-struct','neuromag_meg','-v7');

%%
cd('/group/bgross/work/CIDMEG/Analysis/BrainstormProcess/Brainstorm_db/Pilot2_092025/data/mcwa086_v1/slow')
cd('/group/bgross/work/CIDMEG/Analysis/BrainstormProcess/Brainstorm_db/Pilot2_092025/data/mcwa086_v1/fast')
cd('/group/bgross/work/CIDMEG/Analysis/BrainstormProcess/Brainstorm_db/Pilot2_092025/data/mcwa086_v1/timeout')

%%
% Reference full channel file (has >306 channels) to define which rows are MEG
FULL_CH_FILE = 'channel_vectorview306_acc1.mat';

% Derive the 306 MEG indices (MAG+GRAD) from the full channel file
Lfull = load(FULL_CH_FILE,'Channel');
typesF = string({Lfull.Channel.Type});
iMEGref = find(typesF == "MEG MAG" | typesF == "MEG GRAD");
assert(numel(iMEGref)==306, 'Expected 306 MEG channels in %s, found %d.', FULL_CH_FILE, numel(iMEGref));

files = dir('data_*.mat');
for k = 1:numel(files)
    fn = files(k).name;
    S  = load(fn);                         % F, ChannelFlag, ChannelFile, ...

    nF = size(S.F,1);

    if nF > 306
        % Ensure ChannelFlag length matches before cropping
        cf = S.ChannelFlag(:);
        if numel(cf) < nF, cf(end+1:nF,1) = 1; end
        cf = cf(1:nF);

        % Crop to the 306 MEG rows (keeps original sensor order)
        S.F           = S.F(iMEGref, :);
        S.ChannelFlag = cf(iMEGref);

        if isfield(S,'History')
            S.History(end+1,1:3) = {datestr(now),'crop_to_306','Cropped F/flags to 306 MEG'};
        end

        save(fn,'-struct','S','-v7');      % OVERWRITE
        fprintf('%s: CROPPED %d -> 306 rows\n', fn, nF);

    elseif nF == 306
        % Normalize flags to exactly 306
        cf = S.ChannelFlag(:);
        if numel(cf) < 306, cf(end+1:306,1) = 1; end
        if numel(cf) > 306, cf = cf(1:306); end
        S.ChannelFlag = cf;

        save(fn,'-struct','S','-v7');      % OVERWRITE
        fprintf('%s: already 306; flags normalized\n', fn);

    else
        fprintf('%s: has %d rows (<306), manual check\n', fn, nF);
    end
end

%%
cd('/group/bgross/work/CIDMEG/Analysis/BrainstormProcess/Brainstorm_db/Pilot2_092025/data/mcwa086_v1/Resp2')
cd('/group/bgross/work/CIDMEG/Analysis/BrainstormProcess/Brainstorm_db/Pilot2_092025/data/mcwa086_v1/Resp1')

FULL_CH_FILE = 'channel_vectorview306_acc1.mat';

% Derive the 306 MEG indices (MAG+GRAD) from the full channel file
Lfull = load(FULL_CH_FILE,'Channel');
typesF = string({Lfull.Channel.Type});
iMEGref = find(typesF == "MEG MAG" | typesF == "MEG GRAD");
assert(numel(iMEGref)==306, 'Expected 306 MEG channels in %s, found %d.', FULL_CH_FILE, numel(iMEGref));

files = dir('mcwa*.mat');
for k = 1:numel(files)
    fn = files(k).name;
    S  = load(fn);                         % F, ChannelFlag, ChannelFile, ...

    nF = size(S.F,1);

    if nF > 306
        % Ensure ChannelFlag length matches before cropping
        cf = S.ChannelFlag(:);
        if numel(cf) < nF, cf(end+1:nF,1) = 1; end
        cf = cf(1:nF);

        % Crop to the 306 MEG rows (keeps original sensor order)
        S.F           = S.F(iMEGref, :);
        S.ChannelFlag = cf(iMEGref);

        if isfield(S,'History')
            S.History(end+1,1:3) = {datestr(now),'crop_to_306','Cropped F/flags to 306 MEG'};
        end

        save(fn,'-struct','S','-v7');      % OVERWRITE
        fprintf('%s: CROPPED %d -> 306 rows\n', fn, nF);

    elseif nF == 306
        % Normalize flags to exactly 306
        cf = S.ChannelFlag(:);
        if numel(cf) < 306, cf(end+1:306,1) = 1; end
        if numel(cf) > 306, cf = cf(1:306); end
        S.ChannelFlag = cf;

        save(fn,'-struct','S','-v7');      % OVERWRITE
        fprintf('%s: already 306; flags normalized\n', fn);

    else
        fprintf('%s: has %d rows (<306), manual check\n', fn, nF);
    end
end
