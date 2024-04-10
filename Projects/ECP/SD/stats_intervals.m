PostStim =   [7.21644966e-16, 2];
tlength = 0.3;
Overlap = 0.5;

Overlap1 = 1 - Overlap;
w1 = PostStim(1);
l = tlength;
ov = l * Overlap1;
j = 1;
wi = [];

while w1 + l - ov <= PostStim(2)
    wi(j, :) = [w1, w1 + l];
    j = j + 1;
    w1 = w1 + ov;
end

disp(wi)

%%
% Process: t-test zero [-500ms,1950ms]          H0:(X=0), H1:(X<>0)
bst_process('CallProcess', 'process_test_parametric1', sub_demog_data.sFiles_patn, [], ...
    'timewindow',    [wi(1,1), wi(end,end)], ...
    'scoutsel',      {}, ...
    'scoutfunc',     1, ...  % Mean
    'isnorm',        0, ...
    'avgtime',       0, ...
    'Comment',       '', ...
    'Comment', ['stats_patn_', num2str(wi(1,1)), '_', num2str(wi(end,end))], ...
    'test_type',     'ttest_onesample', ...  % One-sample Student's t-test    X~N(m,s)t = mean(X) ./ std(X) .* sqrt(n)      df=n-1
    'tail',          'one+');  % Two-tailed

%%
% Process: Extract time: [0.000s,2.100s]
% sFiles = bst_process('CallProcess', 'process_extract_time', sub_demog_data.sFiles_ctrl, [], ...
%     'timewindow', [wi(1,1), wi(2,1)-0.01], ...
%     'overwrite',  0);
% sFilesB = db_template('importfile');

for k=1:length(wi)-1
    %     k = 2;
    
    % Process: MEAN: [0.000s,2.100s], abs
    sFiles = bst_process('CallProcess', 'process_average_time', sub_demog_data.sFiles_patn, [], ...
        'timewindow', [wi(k,1), wi(k+1,1)-0.01], ...
        'avg_func',   'mean', ...  % Arithmetic average:  mean(x)
        'overwrite',  0, ...
        'source_abs', 1);
    
    
    % Process: t-test zero [-500ms,1950ms]          H0:(X=0), H1:(X<>0)
    bst_process('CallProcess', 'process_test_parametric1', sFiles, [], ...
        'timewindow',    [wi(k,1), wi(k+1,1)-0.01], ...
        'scoutsel',      {}, ...
        'scoutfunc',     1, ...  % Mean
        'isnorm',        0, ...
        'avgtime',       0, ...
        'Comment',       '', ...
        'Comment', ['stats_patn_', num2str(wi(k,1)), '_', num2str(wi(k+1,1)-0.01)], ...
        'test_type',     'ttest_onesample', ...  % One-sample Student's t-test    X~N(m,s)t = mean(X) ./ std(X) .* sqrt(n)      df=n-1
        'tail',          'one+');  % Two-tailed
    
    %     sFiles1 = [];
    %     for i=1:length(sFiles)
    %         sFiles1{i} = sFiles(i).FileName;
    %     end
    %
    %     % Process: FT t-test unequal fdr [1000ms]          H0:(A=B), H1:(A<>B)
    %     bst_process('CallProcess', 'process_ft_sourcestatistics_VY', sFiles1, sFilesB, ...
    %         'timewindow',     [1, 1], ...
    %         'scoutsel',       {}, ...
    %         'scoutfunc',      1, ...  % Mean
    %         'isabs',          0, ...
    %         'avgtime',        0, ...
    %         'randomizations', 10000, ...
    %         'statistictype',  1, ...  % Independent t-test
    %         'tail',           'one+', ...  % Two-tailed
    %         'correctiontype', 4, ...  % fdr
    %         'minnbchan',      0, ...
    %         'Comment', ['permstats_patn_', num2str(wi(k,1)), '_', num2str(wi(k+1,1)-0.01)], ...
    %         'clusteralpha',   0.05);
    
    dd = rdir('./*abs_mean.mat');
    for i = 1:length(dd)
        disp(dd(i).name)
        pause(0.1),
        delete(dd(i).name)
    end
    
    pause,
    
end