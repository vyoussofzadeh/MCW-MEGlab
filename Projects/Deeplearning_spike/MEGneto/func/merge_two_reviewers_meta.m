function M = merge_two_reviewers_meta(root, r1, r2, Fs)
% root = '/MEG_data/AHW_SpikeAnalysis/Processed/Spike-NoSpike_Mat_files'
% r1,r2 e.g. 'Adi','Josh'

d1 = fullfile(root, r1, 'meta_files');
d2 = fullfile(root, r2, 'meta_files');

f1 = dir(fullfile(d1, ['*_' r1 '_meta_data.txt']));
f2 = dir(fullfile(d2, ['*_' r2 '_meta_data.txt']));

% base name = everything before _<RATER>_meta_data.txt
base1 = erase({f1.name}, ['_' r1 '_meta_data.txt']);
base2 = erase({f2.name}, ['_' r2 '_meta_data.txt']);

common = intersect(base1, base2);

M = struct();
M.common_bases = common;

M.byBase = containers.Map();

for i = 1:numel(common)
    base = common{i};

    meta1_fn = fullfile(d1, [base '_' r1 '_meta_data.txt']);
    meta2_fn = fullfile(d2, [base '_' r2 '_meta_data.txt']);

    meta1 = read_meta_txt(meta1_fn, Fs, false);
    meta2 = read_meta_txt(meta2_fn, Fs, false);

    T1 = meta1.table;  T2 = meta2.table;

    % rename agree columns to keep both
    T1.Properties.VariableNames{'agree'} = ['agree_' r1];
    T2.Properties.VariableNames{'agree'} = ['agree_' r2];

    % merge by candidate identity
    Tm = outerjoin(T1, T2, ...
        'Keys', {'sample','code'}, ...
        'MergeKeys', true, ...
        'Type', 'full');

    % Add useful derived columns
    a1 = Tm.(['agree_' r1]);
    a2 = Tm.(['agree_' r2]);

    bothReviewed = ~isnan(a1) & ~isnan(a2);
    Tm.bothReviewed = bothReviewed;

    Tm.bothAgree    = bothReviewed & (a1==1) & (a2==1);
    Tm.bothReject   = bothReviewed & (a1==0) & (a2==0);
    Tm.disagree     = bothReviewed & (a1~=a2);

    % agreement stats (on overlap only)
    nOverlap = sum(bothReviewed);
    if nOverlap > 0
        Tm_agree_rate = (sum(Tm.bothAgree) + sum(Tm.bothReject)) / nOverlap;
    else
        Tm_agree_rate = NaN;
    end

    M.byBase(base) = struct('T',Tm, 'agreeRate',Tm_agree_rate, ...
                            'meta1_fn',meta1_fn, 'meta2_fn',meta2_fn);
end
end
