%% Analysis of Behavioral Data from MEG Studies

% This script processes behavioral data from MEG experiments conducted at
% the MEG Laboratory, Medical College of Wisconsin & Frodtert Hospital.
% It computes and plots the accuracies for recognizing animals found in the
% US and for matching symbols in the falsefont condition, based on subject
% responses across different runs.

% Author: Vahab Youssof Zadeh
% Date: 5-17-2024
% Version: 1.0
% Usage: This script is intended for batch processing of text files containing
% behavioral data from specific MEG studies. It automatically reads data,
% computes mean accuracies for animal recognition and falsefont matching tasks,
% and generates plots of these accuracies both by subject and across all subjects.

%%
clc; clear; close all

% Base path where the data files are located
basePath = '/group/jbinder/ECP/MEG/MEG_behav/SD/';
savePath = '/group/jbinder/ECP/MEG/MEG_behav/Script';

% Pattern to match files of interest
filePattern = 'SDMatch_eventdata_*.txt';

% Get a list of all files that match the pattern
fileList = dir(fullfile(basePath, filePattern));

% Define the list of common animals
animalsInUS = {'muskrat', 'owl', 'pike', 'mink', 'lamb', 'crow', 'cod', 'bass', 'kitten', 'pelican', 'bluejay', 'robin', 'trout', 'chipmunk', 'perch', 'duck', 'sparrow', 'cardinal', 'cat', 'grasshopper', 'elk', 'goldfish', 'ram', 'turtle', 'salamander', 'parrakeet', 'cockroach', 'horse', 'chicken', 'deer', 'catfish', 'hawk', 'guineapig', 'termite', 'butterfly', 'fox', 'clam', 'caterpillar', 'sow', 'hog', 'bear', 'coyote', 'moose', 'mouse', 'skunk', 'shrimp', 'calf', 'parrot', 'bobcat', 'beaver', 'gnat', 'dolphin', 'rat', 'cow', 'wasp', 'porcupine', 'goose', 'goat', 'spider', 'hamster', 'flounder', 'groundhog', 'dragonfly', 'pig', 'salmon', 'mule', 'lizard', 'carp', 'finch', 'swan', 'ant', 'dog', 'badger', 'sunfish', 'sheep', 'cricket', 'opposum', 'beetle', 'tuna', 'steer', 'chickadee', 'lobster', 'oyster', 'flamingo', 'rabbit', 'moth', 'gerbil', 'boar', 'shrew', 'squirrel', 'pigeon', 'pony', 'turkey', 'cougar', 'stork'};
usedByPeople = {'pike', 'dromedary', 'mink', 'lamb', 'cod', 'bass', 'kitten', 'trout', 'perch', 'duck', 'llama', 'cat', 'goldfish', 'turtle', 'parrakeet', 'horse', 'chicken', 'deer', 'catfish', 'sardine', 'guineapig', 'fox', 'clam', 'sow', 'hog', 'shrimp', 'calf', 'parrot', 'beaver', 'squid', 'cow', 'goose', 'goat', 'hamster', 'flounder', 'pig', 'salmon', 'mule', 'carp', 'dog', 'sunfish', 'sheep', 'ostrich', 'tuna', 'steer', 'lobster', 'camel', 'oyster', 'rabbit', 'gerbil', 'pony', 'turkey'};
commonAnimals = intersect(animalsInUS, usedByPeople);

% Initialize a list to store the extracted information
allSubjectData = [];
acc = struct('animalAccuracy', [], 'falsefontAccuracy', [], 'filename', [], 'subject', [], 'run', []);
accuracyResults = table([], [], [], [], [], [], [], [], 'VariableNames', {'Subject', 'Run', 'Animal_ACC', 'Falsefont_ACC', 'numFontACC', 'numAnimalACC','numAnimalTrials', 'numFalsefontTrials'});

% Loop through each file in the file list
for iFile = 1:length(fileList)
    filename = fileList(iFile).name;
    % Extract subject and run from the filename
    tokens = regexp(filename, 'SDMatch_eventdata_(\d+)_(\d+).txt', 'tokens');
    subject = str2double(tokens{1}{1});
    run = str2double(tokens{1}{2});
    
    % Full path to the file
    fullFilename = fullfile(basePath, filename);
    
    % Read the data file
    data = readtable(fullFilename);
    
    % Initialize a list to store the extracted information for this file
    extractedData = [];
    isExperimentStarted = false;
    
    % Loop through each row in the data table to extract information
    for i = 1:height(data)
        message = data.Var3{i};  % Get the current message
        
        % Check if the experiment has started
        if contains(message, 'Keypress')
            isExperimentStarted = true;
        end
        
        % Process messages only after the experiment has started
        if isExperimentStarted
            % Initialize variables for the current trial
            currentCondition = "";
            currentResp = "";
            currentFontCresp = "";
            currentFontACC = "";
            currentFontType = "";
            currentAnimalName = "";
            currentAnimalACC = "";
            currentACC_commonAnimals = "";
            
            % Extract condition
            conditionTokens = regexp(message, '"condition":\s*"?(\w+)"?', 'tokens');
            if ~isempty(conditionTokens) && ~isempty(conditionTokens{1})
                currentCondition = conditionTokens{1}{1};
            end
            
            % Extract response
            respTokens = regexp(message, '"keyboard_resp":\s*\[\[(\d+),\s*(\d+)\]\]', 'tokens');
            if ~isempty(respTokens) && ~isempty(respTokens{1})
                currentResp = "1";
            end
            
            % Extract Font_Cresp
            fontCrespTokens = regexp(message, '"Font_Cresp":\s*"?(\w+)"?', 'tokens');
            if ~isempty(fontCrespTokens) && ~isempty(fontCrespTokens{1})
                currentFontCresp = fontCrespTokens{1}{1};
            end
            
            % Extract Font_ACC
            fontACCTokens = regexp(message, '"Font_ACC":\s*"?(\w+)"?', 'tokens');
            if ~isempty(fontACCTokens) && ~isempty(fontACCTokens{1})
                currentFontACC = fontACCTokens{1}{1};
            end
            
            % Extract Font_type (text of created TextStim objects)
            fontTypeTokens = regexp(message, '"stim":\s*"(\w+)"', 'tokens');
            if ~isempty(fontTypeTokens) && ~isempty(fontTypeTokens{1})
                currentFontType = fontTypeTokens{1}{1};
            end
            
            % Extract Animal_name (name of created TextStim objects)
            animalNameTokens = regexp(message, '"stim":\s*"(\w+)"', 'tokens');
            if ~isempty(animalNameTokens) && ~isempty(animalNameTokens{1})
                currentAnimalName = animalNameTokens{1}{1};
            end
            
            % Extract Animal_ACC
            animalACCTokens = regexp(message, '"Animal_ACC":\s*"?(\w+)"?', 'tokens');
            if ~isempty(animalACCTokens) && ~isempty(animalACCTokens{1})
                currentAnimalACC = animalACCTokens{1}{1};
            end
            
            % Set Animal_name to blank if Condition is falsefont
            if currentCondition == "falsefont"
                currentAnimalName = "";
            end
            
            % Set Font_type to blank if Condition is animal
            if currentCondition == "animal"
                currentFontType = "";
            end
            
            % Check for non-empty par_resp and set currentResp to 1 if condition is falsefont
            if currentCondition == "falsefont"
                parRespTokens = regexp(message, '"par_resp":\s*\[([^\]]+)\]', 'tokens');
                if ~isempty(parRespTokens)
                    currentResp = "1";
                end
                
                % Check for num_targets and set currentFontCresp to 1 if num_targets is 2
                numTargetsTokens = regexp(message, '"num_targets":\s*(2)', 'tokens');
                if ~isempty(numTargetsTokens)
                    currentFontCresp = "1";
                end
            end
            
            % Check for non-empty par_resp and set currentResp to 1 if condition is animal
            if currentCondition == "animal"
                parRespTokens = regexp(message, '"par_resp":\s*\[([^\]]+)\]', 'tokens');
                if ~isempty(parRespTokens)
                    currentResp = "1";
                end
            end
            
            % Set Font_ACC to 1 if both resp and Font_Cresp are 1
            if currentResp == "1" && currentFontCresp == "1"
                currentFontACC = "1";
            end
            
            % Check if Animal_name is in commonAnimals and set ACC_commonAnimals to 1
            if any(strcmp(currentAnimalName, commonAnimals))
                currentACC_commonAnimals = "1";
            end
            
            % Set Animal_ACC to 1 if both ACC_commonAnimals and resp are 1
            if currentACC_commonAnimals == "1" && currentResp == "1"
                currentAnimalACC = "1";
            end
            
            % Check if any field is non-empty
            if currentCondition ~= "" || currentResp ~= "" || currentFontCresp ~= "" || ...
                    currentFontACC ~= "" || currentFontType ~= "" || currentAnimalName ~= "" || ...
                    currentAnimalACC ~= "" || currentACC_commonAnimals ~= ""
                % Append the current trial data to the extractedData list
                extractedData = [extractedData; {subject, run, i, currentCondition, currentResp, currentFontCresp, ...
                    currentFontACC, currentFontType, currentAnimalName, currentAnimalACC, currentACC_commonAnimals}];
            end
        end
    end
    
    % Create the table with the extracted information
    trialData = cell2table(extractedData, ...
        'VariableNames', {'Subject', 'Run', 'Trial', 'Condition', 'resp', 'Font_Cresp', 'Font_ACC', 'Font_type', 'Animal_name', 'Animal_ACC', 'ACC_commonAnimals'});
    
    numFalsefontTrials_all = sum(trialData.resp == "" & trialData.Condition == "falsefont");
    numAnimalTrials_all = sum(trialData.ACC_commonAnimals == "" & trialData.Condition == "animal");
    
    numFalsefontTrials = sum(trialData.Font_Cresp == "1" & trialData.Condition == "falsefont");
    numAnimalTrials = sum(trialData.ACC_commonAnimals == "1" & trialData.Condition == "animal");
    
    numFontACC = sum(trialData.Font_ACC == "1" & trialData.Condition == "falsefont");
    numAnimalACC = sum(trialData.Animal_ACC == "1" & trialData.Condition == "animal");
    
    falsefontAccuracy = (numFontACC / numFalsefontTrials) * 100;
    animalAccuracy = (numAnimalACC / numAnimalTrials) * 100;
    
    fprintf('Subject %d Run %d - Falsefont Accuracy: %.2f%%\n', subject, run, falsefontAccuracy);
    fprintf('Subject %d Run %d - Animal Accuracy: %.2f%%\n', subject, run, animalAccuracy);
    
    acc.animalAccuracy(iFile) = animalAccuracy;
    acc.falsefontAccuracy(iFile) = falsefontAccuracy;
    acc.filename{iFile} = filename;
    acc.subject(iFile) = subject;
    acc.run(iFile) = run;
    
    % Append results to the accuracyResults table
    newResult = {subject, run, animalAccuracy, falsefontAccuracy, numFontACC, numAnimalACC, numAnimalTrials, numFalsefontTrials};
    accuracyResults = [accuracyResults; newResult];    
    fprintf('Processed file: %s\n', fullFilename);
end

% Display the summary table
disp(accuracyResults);

%% Plotting 
meanAnimalAcc = nanmean(accuracyResults.Animal_ACC);
meanFalsefontAcc = nanmean(accuracyResults.Falsefont_ACC);
fprintf('Mean Accuracy for Animal Condition: %.2f%%\n', meanAnimalAcc);
fprintf('Mean Accuracy for Falsefont Condition: %.2f%%\n', meanFalsefontAcc);

meanAnimaltargets = nanmean(accuracyResults.numAnimalTrials)/64;
meanSymboltargets = nanmean(accuracyResults.numFalsefontTrials)/64;
fprintf('Mean targets for Animal Condition: %.2f%%\n', meanAnimaltargets);
fprintf('Mean targets for Falsefont Condition: %.2f%%\n', meanSymboltargets);

% Calculate mean accuracy per subject for each condition
meanAccBySubject_Animal = groupsummary(accuracyResults, 'Subject', 'mean', 'Animal_ACC');
meanAccBySubject_Falsefont = groupsummary(accuracyResults, 'Subject', 'mean', 'Falsefont_ACC');

L = length(meanAccBySubject_Animal.Subject);

% Plot for Animal Condition
figure;
subplot 211
bar(meanAccBySubject_Animal.mean_Animal_ACC)
title({'Anim Acc'});
xlabel('Subject');
ylabel('Mean Accuracy');
ylim([0 100]); % Assuming accuracy is a proportion
set(gca,'Xtick', 1:L,'XtickLabel',num2str(meanAccBySubject_Animal.Subject));
xtickangle(90); % Angle the subject labels for better readability

% Plot for Falsefont Condition
% figure;
subplot 212
bar(meanAccBySubject_Falsefont.mean_Falsefont_ACC);
title({'FalseFont Acc'});
xlabel('Subject');
ylabel('Mean Accuracy');
ylim([0 100]); % Assuming accuracy is a proportion
set(gca,'Xtick', 1:L,'XtickLabel',num2str(meanAccBySubject_Falsefont.Subject));
xtickangle(90); % Angle the subject labels for better readability

%% Mean and std summary
% Overall mean and standard deviation for Animal condition
overallMean_Animal = nanmean(meanAccBySubject_Animal.mean_Animal_ACC);
stdError_Animal = nanstd(meanAccBySubject_Animal.mean_Animal_ACC) / sqrt(height(meanAccBySubject_Animal));

% Overall mean and standard deviation for Falsefont condition
overallMean_Falsefont = nanmean(meanAccBySubject_Falsefont.mean_Falsefont_ACC);
stdError_Falsefont = nanstd(meanAccBySubject_Falsefont.mean_Falsefont_ACC) / sqrt(height(meanAccBySubject_Falsefont));

%% Plot mean and std
% X-coordinates for the two conditions
conditions = [1, 2]; 

% Mean accuracies
means = [overallMean_Animal, overallMean_Falsefont];

% Standard errors
errors = [stdError_Animal, stdError_Falsefont];

% Creating the error bar plot
figure;
errorbar(conditions, means, errors, 'o');
xlim([0.8 2.2]); % Adjust x-axis limits for better visualization
xticks(conditions);
xticklabels({'Animal', 'Falsefont'});
ylabel('Mean Accuracy');
title('Mean Accuracy with Standard Error for Each Condition');
axis square

%% Save as mat file
% Define the file path where you want to save the findings
fileName = 'TaskPerformanceSD.mat'; % Name of the .mat file

% Full path to the .mat file
fullFilePath = fullfile(savePath, fileName);

% Save specific variables to the .mat file
save(fullFilePath, 'accuracyResults', 'meanAnimalAcc', 'meanFalsefontAcc', 'overallMean_Animal', 'overallMean_Falsefont', 'stdError_Animal', 'stdError_Falsefont');


%% Save results as CSV file
csvFileName = 'TaskPerformanceSD.csv'; % Name of the CSV file
csvFilePath = fullfile(savePath, csvFileName);

% Save the accuracy results table as a CSV file
writetable(accuracyResults, csvFilePath, 'Delimiter', ',');