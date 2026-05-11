function S = get_subject_epochs(trials, subj)
% subj: string or char, exact match to trials.groups entries
mask = strcmp(trials.groups, subj);

S.subj      = subj;
S.idx       = find(mask);
S.data      = trials.data(mask);
S.labels    = trials.labels(mask);
% S.refLbl    = trials.refLbl;

% class masks within the subject
S.isSpike   = (S.labels=='Spike');
S.isNoSpike = (S.labels=='NoSpike');

% split by class
S.XSpike    = S.data(S.isSpike);
S.XNoSpike  = S.data(S.isNoSpike);

% counts
S.nSpike    = nnz(S.isSpike);
S.nNoSpike  = nnz(S.isNoSpike);
S.nTotal    = numel(S.data);
end
