function patn_MEGfMRI_neuropsych = ecpfunc_sum_neuropsych(sub_MF_pt)

patn_neuropsych_tle = ecpfunc_read_patn_neuropsych_tle();
% TLESide = patn_neuropsych_tle.TLESide; 
SUBNO = patn_neuropsych_tle.SUBNO;

% clc
SUBNO_pt2 = [];
for i=1:length(sub_MF_pt)
    SUBNO_pt2(i) = str2double(sub_MF_pt{i}(3:end));
end

[~,~,IB_tle] = intersect(SUBNO_pt2, SUBNO);
size(IB_tle);
patn_MEGfMRI_neuropsych = patn_neuropsych_tle(IB_tle,:);

% summaryStats = @(x) table(mean(x,'omitnan'), median(x,'omitnan'), std(x,'omitnan'));
summaryStats = @(x) table(mean(x,'omitnan'), median(x,'omitnan'), std(x,'omitnan'), ...
                          'VariableNames', {'Mean', 'Median', 'Std'});

stats = varfun(summaryStats, patn_MEGfMRI_neuropsych, 'InputVariables', ...
    {'Age', 'EHQ', 'EducationYRS', 'WASI_BlckR_ZScore', 'WASI_VocR_ZScore'}, ...
    'OutputFormat', 'table');
disp(stats);

variablesToCorrelate = {'Age', 'EHQ', 'EducationYRS', 'WASI_BlckR_ZScore', 'WASI_VocR_ZScore'};
correlationMatrix = corr(patn_MEGfMRI_neuropsych{:, variablesToCorrelate}, 'Rows', 'complete');
disp(array2table(correlationMatrix, 'VariableNames', variablesToCorrelate, 'RowNames', variablesToCorrelate));

lm = fitlm(patn_MEGfMRI_neuropsych, 'Age~EHQ');
disp(lm);

% Gender summary
genderSummary = groupsummary(patn_MEGfMRI_neuropsych, 'Sex', 'IncludeMissingGroups', true);
disp(genderSummary);

% TLE side summary
TLESummary = groupsummary(patn_MEGfMRI_neuropsych, 'TLESide', 'IncludeMissingGroups', true);
disp(TLESummary);

% Age statistics
meanAge = mean(patn_MEGfMRI_neuropsych.Age, 'omitnan');
stdAge = std(patn_MEGfMRI_neuropsych.Age, 'omitnan');
disp(['Mean Age: ', num2str(meanAge)]);
disp(['Standard Deviation of Age: ', num2str(stdAge)]);

minAge = min(patn_MEGfMRI_neuropsych.Age);
maxAge = max(patn_MEGfMRI_neuropsych.Age);
ageRange = maxAge - minAge;
disp(['Minimum Age: ', num2str(minAge)]);
disp(['Maximum Age: ', num2str(maxAge)]);
disp(['Age Range: ', num2str(ageRange)]);

% Handedness summary
handednessSummary = groupsummary(patn_MEGfMRI_neuropsych, 'DomntHand', 'IncludeMissingGroups', true);
disp(handednessSummary);

