%% Analysis of Behavioral Data from MEG Studies

% This script processes behavioral data from MEG experiments conducted at
% the MEG Laboratory, Medical College of Wisconsin & Frodtert Hospital.
% It computes and plots the accuracies for recognizing animals found in the
% US and for matching symbols in the falsefont condition, based on subject
% responses across different runs.

% Author: Vahab Youssof Zadeh
% Date: 4-10-2024
% Version: 1.0
% Usage: This script is intended for batch processing of text files containing
% behavioral data from specific MEG studies. It automatically reads data,
% computes mean accuracies for animal recognition and falsefont matching tasks,
% and generates plots of these accuracies both by subject and across all subjects.

clc; clear; close all

% Base path where the data files are located
basePath = '/group/jbinder/ECP/MEG/MEG_behav/SD_score/';
savePath = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data';

% Pattern to match files of interest
filePattern = 'SDMatch_score_*.txt';

% Get a list of all files that match the pattern
fileList = dir(fullfile(basePath, filePattern));

sd_correct.found_in_US = {'muskrat', 'owl', 'pike', 'mink', 'lamb', 'crow', 'cod', 'bass', 'kitten', 'pelican', 'bluejay', 'robin', 'trout', 'chipmunk', 'perch', 'duck', 'sparrow', 'cardinal', 'cat', 'grasshopper', 'elk', 'goldfish', 'ram', 'turtle', 'salamander', 'parrakeet', 'cockroach', 'horse', 'chicken', 'deer', 'catfish', 'hawk', 'guineapig', 'termite', 'butterfly', 'fox', 'clam', 'caterpillar', 'sow', 'hog', 'bear', 'coyote', 'moose', 'mouse', 'skunk', 'shrimp', 'calf', 'parrot', 'bobcat', 'beaver', 'gnat', 'dolphin', 'rat', 'cow', 'wasp', 'porcupine', 'goose', 'goat', 'spider', 'hamster', 'flounder', 'groundhog', 'dragonfly', 'pig', 'salmon', 'mule', 'lizard', 'carp', 'finch', 'swan', 'ant', 'dog', 'badger', 'sunfish', 'sheep', 'cricket', 'opposum', 'beetle', 'tuna', 'steer', 'chickadee', 'lobster', 'oyster', 'flamingo', 'rabbit', 'moth', 'gerbil', 'boar', 'shrew', 'squirrel', 'pigeon', 'pony', 'turkey', 'cougar', 'stork'};
sd_correct.used_by_people = {'pike', 'dromedary', 'mink', 'lamb', 'cod', 'bass', 'kitten', 'trout', 'perch', 'duck', 'llama', 'cat', 'goldfish', 'turtle', 'parrakeet', 'horse', 'chicken', 'deer', 'catfish', 'sardine', 'guineapig', 'fox', 'clam', 'sow', 'hog', 'shrimp', 'calf', 'parrot', 'beaver', 'squid', 'cow', 'goose', 'goat', 'hamster', 'flounder', 'pig', 'salmon', 'mule', 'carp', 'dog', 'sunfish', 'sheep', 'ostrich', 'tuna', 'steer', 'lobster', 'camel', 'oyster', 'rabbit', 'gerbil', 'pony', 'turkey'};
sd_correct.four_legs = {'chimpanzee', 'gorilla', 'ape', 'muskrat', 'dromedary', 'mink', 'lamb', 'kitten', 'lion', 'chipmunk', 'reindeer', 'llama', 'cat', 'elk', 'crocodile', 'tortoise', 'ram', 'turtle', 'salamander', 'horse', 'deer', 'jackal', 'guineapig', 'tiger', 'leopard', 'fox', 'elephant', 'sow', 'hog', 'bear', 'coyote', 'moose', 'mouse', 'skunk', 'rhinocerous', 'jaguar', 'calf', 'bobcat', 'hyena', 'beaver', 'rat', 'cow', 'porcupine', 'iguana', 'gazelle', 'zebra', 'goat', 'anteater', 'hamster', 'groundhog', 'monkey', 'pig', 'mule', 'lizard', 'giraffe', 'cheetah', 'panther', 'dog', 'badger', 'sloth', 'sheep', 'opposum', 'steer', 'camel', 'rabbit', 'gerbil', 'boar', 'shrew', 'squirrel', 'pony', 'cougar'};


% Initialize a table to store accuracy results
accuracyResults = table([], [], [], [], 'VariableNames', {'Subject', 'Run', 'Animal_ACC', 'Falsefont_ACC'});

% Loop through each file found
for iFile = 1:length(fileList)
    filename = fileList(iFile).name;
    % Extract subject and run from the filename
    tokens = regexp(filename, 'SDMatch_score_(\d+)_(\d+).txt', 'tokens');
    subject = str2double(tokens{1}{1});
    run = str2double(tokens{1}{2});
    
    % Full path to the file
    fullFilename = fullfile(basePath, filename);
    
    % Read the table from the file
    opts = detectImportOptions(fullFilename, 'FileType', 'text');
    opts.Delimiter = '\t';
    data = readtable(fullFilename, opts);
    
    % Rename columns for clarity
    data.Properties.VariableNames = {'Trial', 'Condition', 'Response', 'Font_Cresp', 'Font_ACC', 'Font_type', 'Animal_name', 'ExtraVar1'};
    
    % Compute accuracy for animal and falsefont conditions
    % Filter for 'animal' condition trials
    animalTrials = data(strcmp(data.Condition, 'animal'), :);
    
    % Initialize count of correct responses for animals found in the US
    correctAnimalCount = 0;
    
    % Loop through each 'animal' trial
    for i = 1:height(animalTrials)
        animalName = animalTrials.Animal_name{i};
        response = animalTrials.Response(i); 
        % Check if the response is numeric and convert to string if necessary
        if isnumeric(response)
            response = num2str(response); % Converts numeric to character array
        end
        
        % Check if the animal is in the list of animals found in the US
        isInUS = ismember(animalName, sd_correct.found_in_US);
        
        % Assuming "1" means the participant identified the animal correctly
        if isInUS && strcmp(response, '1')
            correctAnimalCount = correctAnimalCount + 1;
        elseif ~isInUS && ~strcmp(response, '1') % If specific response for not found in the US is considered
            correctAnimalCount = correctAnimalCount + 1;
        end
    end
    
    % Calculate accuracy for the animal condition
    animalAcc = correctAnimalCount / height(animalTrials);
    
    % Filter for 'falsefont' condition trials
    falsefontTrials = data(strcmp(data.Condition, 'falsefont'), :);
    
    % Initialize count of correct responses for symbol matching
    correctFalsefontCount = 0;
    
    % Loop through each 'falsefont' trial
    for i = 1:height(falsefontTrials)
        response = falsefontTrials.Response(i); 
        
        % Check if the response is numeric and convert to string if necessary
        if isnumeric(response)
            response = num2str(response); % Converts numeric to character array
        end
%         disp(response)
        fontAcc = falsefontTrials.Font_ACC(i);
        % Check if the participant's response ("1" for a match) is correct based on Font_ACC
        if strcmp(response, '1') && fontAcc == 1
            correctFalsefontCount = correctFalsefontCount + 1;
        elseif (~strcmp(response, '1') && fontAcc == 0) % Consider non-"1" responses as correct for non-matches
            correctFalsefontCount = correctFalsefontCount + 1;
        end
    end
    
    % Calculate accuracy for the falsefont condition
    falsefontAcc = correctFalsefontCount / height(falsefontTrials);
    
    
    % Append results to the accuracyResults table
    newResult = {subject, run, animalAcc, falsefontAcc};
    accuracyResults = [accuracyResults; newResult];
    
    fprintf('Processed file: %s\n', fullFilename);
end

disp(accuracyResults)
% Summarize accuracy values
meanAnimalAcc = mean(accuracyResults.Animal_ACC);
meanFalsefontAcc = mean(accuracyResults.Falsefont_ACC);
fprintf('Mean Accuracy for Animal Condition: %.2f%%\n', meanAnimalAcc * 100);
fprintf('Mean Accuracy for Falsefont Condition: %.2f%%\n', meanFalsefontAcc * 100);


% Calculate mean accuracy per subject for each condition
meanAccBySubject_Animal = groupsummary(accuracyResults, 'Subject', 'mean', 'Animal_ACC');
meanAccBySubject_Falsefont = groupsummary(accuracyResults, 'Subject', 'mean', 'Falsefont_ACC');

L = length(meanAccBySubject_Animal.Subject);

% Plot for Animal Condition
figure;
% bar(meanAccBySubject_Animal.Subject);
bar(meanAccBySubject_Animal.mean_Animal_ACC)
title('Mean Accuracy for Animal Condition by Subject');
xlabel('Subject');
ylabel('Mean Accuracy');
ylim([0 1]); % Assuming accuracy is a proportion
set(gca,'Xtick', 1:L,'XtickLabel',num2str(meanAccBySubject_Animal.Subject));
xtickangle(90); % Angle the subject labels for better readability

% Plot for Falsefont Condition
figure;
bar(meanAccBySubject_Falsefont.mean_Falsefont_ACC);
title('Mean Accuracy for Falsefont Condition by Subject');
xlabel('Subject');
ylabel('Mean Accuracy');
ylim([0 1]); % Assuming accuracy is a proportion
set(gca,'Xtick', 1:L,'XtickLabel',num2str(meanAccBySubject_Falsefont.Subject));
xtickangle(90); % Angle the subject labels for better readability

%%
% Overall mean and standard deviation for Animal condition
overallMean_Animal = mean(meanAccBySubject_Animal.mean_Animal_ACC);
stdError_Animal = std(meanAccBySubject_Animal.mean_Animal_ACC) / sqrt(height(meanAccBySubject_Animal));

% Overall mean and standard deviation for Falsefont condition
overallMean_Falsefont = mean(meanAccBySubject_Falsefont.mean_Falsefont_ACC);
stdError_Falsefont = std(meanAccBySubject_Falsefont.mean_Falsefont_ACC) / sqrt(height(meanAccBySubject_Falsefont));

%%
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

%%
% Define the file path where you want to save the findings
fileName = 'SD_taskperformace.mat'; % Name of the .mat file

% Full path to the .mat file
fullFilePath = fullfile(savePath, fileName);

% Save specific variables to the .mat file
save(fullFilePath, 'accuracyResults', 'meanAnimalAcc', 'meanFalsefontAcc', 'overallMean_Animal', 'overallMean_Falsefont', 'stdError_Animal', 'stdError_Falsefont');


