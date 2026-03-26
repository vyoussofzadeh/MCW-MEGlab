function plot_two_rater_data_categories_fromX(Xs, Tm, r1, r2, Fs)
% Xs: N×C×T, Tm has category flags: bothAgree/bothReject/disagree

    if nargin < 6 || isempty(Fs) || isnan(Fs) || Fs<=0
        Fs = 1; % avoid crash; x-axis becomes samples
    end
    t = (0:size(Xs,3)-1)/Fs;

    bothAgree  = Tm.bothAgree;
    bothReject = Tm.bothReject;
    disagree   = Tm.disagree;

    mA = squeeze(mean(Xs(bothAgree,:,:),  1, 'omitnan'));
    mR = squeeze(mean(Xs(bothReject,:,:), 1, 'omitnan'));
    mD = squeeze(mean(Xs(disagree,:,:),   1, 'omitnan'));

    figure; imagesc(t, 1:size(Xs,2), mA); axis tight; colorbar;
    title(sprintf('%s vs %s: BOTH AGREE mean', r1, r2)); xlabel('Time (s)'); ylabel('Channel');

    figure; imagesc(t, 1:size(Xs,2), mR); axis tight; colorbar;
    title(sprintf('%s vs %s: BOTH REJECT mean', r1, r2)); xlabel('Time (s)'); ylabel('Channel');

    figure; imagesc(t, 1:size(Xs,2), mD); axis tight; colorbar;
    title(sprintf('%s vs %s: DISAGREE mean', r1, r2)); xlabel('Time (s)'); ylabel('Channel');

    % example traces
    ch = 1; nExamples = 5;
    figure; hold on;
    idxA = find(bothAgree);  idxA = idxA(1:min(nExamples,numel(idxA)));
    idxR = find(bothReject); idxR = idxR(1:min(nExamples,numel(idxR)));
    idxD = find(disagree);   idxD = idxD(1:min(nExamples,numel(idxD)));

    for i = 1:numel(idxA), plot(t, squeeze(Xs(idxA(i),ch,:)), 'g-'); end
    for i = 1:numel(idxR), plot(t, squeeze(Xs(idxR(i),ch,:)), 'r-'); end
    for i = 1:numel(idxD), plot(t, squeeze(Xs(idxD(i),ch,:)), 'k-'); end
    grid on; xlabel('Time (s)'); ylabel('Amplitude');
    title(sprintf('Example traces ch=%d (g both agree, r both reject, k disagree)', ch));
end
