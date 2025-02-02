taskperf_datadir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/';
sub_MF_pt_num = cellfun(@(x) str2double(x(3:end)), sub_MF_pt);
taskPerformanceDataPath = fullfile(taskperf_datadir, 'TaskPerformanceSD.mat');
load(taskPerformanceDataPath);

% Calculating the average performance across runs while handling missing data
uniqueSubjects = unique(accuracyResults.Subject);
avgAnimalAcc = zeros(length(uniqueSubjects), 1);
avgFalsefontAcc = zeros(length(uniqueSubjects), 1);

for i = 1:length(uniqueSubjects)
    subjectData = accuracyResults(accuracyResults.Subject == uniqueSubjects(i), :);
    
    % Filter out 0 or NaN values for Animal_ACC
    validAnimalScores = subjectData.Animal_ACC(subjectData.Animal_ACC > 0 & ~isnan(subjectData.Animal_ACC));
    if ~isempty(validAnimalScores)
        avgAnimalAcc(i) = mean(validAnimalScores);
    else
        avgAnimalAcc(i) = NaN; % Assign NaN if all scores were invalid
    end
    
    % Filter out 0 or NaN values for Falsefont_ACC
    validFalsefontScores = subjectData.Falsefont_ACC(subjectData.Falsefont_ACC > 0 & ~isnan(subjectData.Falsefont_ACC));
    if ~isempty(validFalsefontScores)
        avgFalsefontAcc(i) = mean(validFalsefontScores);
    else
        avgFalsefontAcc(i) = NaN; % Assign NaN if all scores were invalid
    end
end

% Create a table to store the results
avgResults = table(uniqueSubjects, avgAnimalAcc, avgFalsefontAcc, ...
                   'VariableNames', {'SubjectID', 'Animal_ACC', 'Symbol_ACC'});
               
[~, ~, IB_taskperformance] = intersect(sub_MF_pt_num, avgResults.SubjectID);

% Extract updated accuracy results for the subjects of interest
accuracyResults_updt = avgResults(IB_taskperformance,:);

% Display overall statistics
totalmean = mean(avgResults.Symbol_ACC, 'omitnan');
disp(['The total mean of Symbol_ACC is: ', num2str(totalmean)]);
totalstd = std(avgResults.Symbol_ACC, 'omitnan');
disp(['The total std of Symbol_ACC is: ', num2str(totalstd)]);
totalmean = mean(avgResults.Animal_ACC, 'omitnan');
disp(['The total mean of mean_Animal_ACC is: ', num2str(totalmean)]);
totalstd = std(avgResults.Animal_ACC, 'omitnan');
disp(['The total std of mean_Animal_ACC is: ', num2str(totalstd)]);

% Update SubjectID format to include 'EC' prefix
accuracyResults_updt.SubjectID = cellstr(strcat('EC', string(accuracyResults_updt.SubjectID)));