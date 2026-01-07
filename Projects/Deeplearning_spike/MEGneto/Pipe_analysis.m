addpath('/MEG_data/LAB_MEMBERS/Vahab/Github/MCW-MEGlab/MCW_MEGlab_git/Projects/Deeplearning_spike/function')

root = '/MEG_data/AHW_SpikeAnalysis/Processed/Spike-NoSpike_Mat_files';
raters = {'Adi','Josh','Manoj','Pradeep'};

allX = [];
ally = [];
pairs = {};

for r = 1:numel(raters)
    d = fullfile(root, raters{r});
    spk = dir(fullfile(d, '*_spike_events.mat'));
    for k = 1:numel(spk)
        spk_fn = fullfile(spk(k).folder, spk(k).name);

        base = erase(spk(k).name, '_spike_events.mat');
        nos_fn = fullfile(d, [base '_nospike_events.mat']);

        if ~isfile(nos_fn)
            fprintf('Missing nospike pair for: %s\n', spk(k).name);
            continue;
        end

        [Es, metas] = load_epochs_mat(spk_fn);
        [En, metan] = load_epochs_mat(nos_fn);

        % Expect: trials × time × chan  (like your 19×343×21)
        % Convert to: trials × chan × time
        Xs = permute(Es(:,1:306,:), [1 3 2]);
        Xn = permute(En(:,1:306,:), [1 3 2]);

        X = cat(1, Xs, Xn);
        y = [ones(size(Xs,1),1); zeros(size(Xn,1),1)];

        allX = cat(1, allX, X);
        ally = cat(1, ally, y);

        pairs(end+1,:) = {raters{r}, base, spk_fn, nos_fn, size(Xs,1), size(Xn,1), size(Xs,2), size(Xs,3)}; % #ok<AGROW>

        % (Optional) sanity check Fs consistency if present
        if isfield(metas,'Fs') && isfield(metan,'Fs') && metas.Fs ~= metan.Fs
            warning('Fs mismatch for %s (%g vs %g)', base, metas.Fs, metan.Fs);
        end
    end
end

T = cell2table(pairs, 'VariableNames', {'Rater','Base','SpikeFile','NoSpikeFile','nSpike','nNoSpike','nChan','nTime'});
disp(T)

fprintf('Final dataset: X = [%d trials × %d chan × %d time], y = [%d×1]\n', ...
    size(allX,1), size(allX,2), size(allX,3), numel(ally));

%%
Fs = 128; % set if not stored; otherwise use meta.Fs when available
trial = 1;
t = (0:size(allX,3)-1)/Fs;

figure;
plot(t, squeeze(allX(trial,1,:)));
xlabel('Time (s)'); ylabel('Amplitude');
title(sprintf('Trial %d, Channel 1, label=%d', trial, ally(trial)));
grid on;

figure;
imagesc(t, 1:size(allX,2), squeeze(allX(trial,:,:)));
axis tight;
xlabel('Time (s)'); ylabel('Channel');
title('Trial heatmap (channels × time)');
colorbar;
