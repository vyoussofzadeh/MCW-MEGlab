%% ECP Semantic Decision Task Dataset, Medical College of Wisconsin

% Script: BS Process (Laterality Analysis)
% Project: ECP_SD
% Written by: Vahab Youssof Zadeh

clear; clc; close all; warning off;

%% Paths
restoredefaultpath
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/run')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/function')
Run_setpath
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/External/other/')
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/tools/helpful_tools/daviolinplot/daboxplot');

%% 72 patients
sub_MF_pt = {'EC1007'
    'EC1014'
    'EC1020'
    'EC1028'
    'EC1029'
    'EC1030'
    'EC1042'
    'EC1047'
    'EC1048'
    'EC1060'
    'EC1067'
    'EC1068'
    'EC1072'
    'EC1073'
    'EC1074'
    'EC1075'
    'EC1077'
    'EC1078'
    'EC1079'
    'EC1081'
    'EC1082'
    'EC1097'
    'EC1099'
    'EC1100'
    'EC1101'
    'EC1102'
    'EC1103'
    'EC1106'
    'EC1109'
    'EC1113'
    'EC1114'
    'EC1115'
    'EC1116'
    'EC1117'
    'EC1118'
    'EC1121'
    'EC1122'
    'EC1123'
    'EC1125'
    'EC1129'
    'EC1130'
    'EC1133'
    'EC1135'
    'EC1138'
    'EC1139'
    'EC1141'
    'EC1142'
    'EC1144'
    'EC1150'
    'EC1153'
    'EC1154'
    'EC1156'
    'EC1157'
    'EC1158'
    'EC1159'
    'EC1161'
    'EC1162'
    'EC1165'
    'EC1167'
    'EC2029'
    'EC2038'
    'EC2045'
    'EC2052'
    'EC2054'
    'EC2072'
    'EC2074'
    'EC2083'
    'EC2085'
    'EC2090'
    'EC2109'
    'EC2112'
    'EC2114'};

%% Response (Reaction) Time Data
load('/data/MEG/Research/aizadi/process/RT_summary/ResponseTime.mat')
[~,~,IB_reactiontime] = intersect(sub_MF_pt, T.Sub_ID);
T_patn_MEGfMRI = T(IB_reactiontime,:);
meanAnimal = mean(T_patn_MEGfMRI.Animal, 'omitnan');
stdAnimal = std(T_patn_MEGfMRI.Animal, 'omitnan');
meanSymbol = mean(T_patn_MEGfMRI.Symbol, 'omitnan');
stdSymbol = std(T_patn_MEGfMRI.Symbol, 'omitnan');
disp(['Mean of Animal reaction times: ', num2str(meanAnimal)]);
disp(['Standard Deviation of Animal reaction times: ', num2str(stdAnimal)]);
disp(['Mean of Symbol reaction times: ', num2str(meanSymbol)]);
disp(['Standard Deviation of Symbol reaction times: ', num2str(stdSymbol)]);
[correlationCoefficient, p] = corr(T_patn_MEGfMRI.Animal, T_patn_MEGfMRI.Symbol, 'Rows', 'complete');
validPairs = sum(~isnan(T_patn_MEGfMRI.Animal) & ~isnan(T_patn_MEGfMRI.Symbol));
df = validPairs - 2;
disp(['Correlation coefficient: ', num2str(correlationCoefficient)]);
disp(['Degrees of freedom: ', num2str(df)]);

%% Response Time Data (2)
figure;
daboxplot(T_patn_MEGfMRI.Avg, 'groups', ones(1, numel(T_patn_MEGfMRI.Avg)), 'outsymbol', 'kx', 'xtlabels', 'test', 'fill', 0);
xlabel('All Subjects');
ylabel('Response Time (sec)');
set(gca, 'FontSize', 10);
set(gca,'color','none');
hold off;
title({'Response Time'});
set(gca, 'XTick', []);
box off;
set(gcf, 'Position', [1000   100   300   300]);

%% Task Performance
taskperf_datadir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/';
sub_MF_pt_num = cellfun(@(x) str2double(x(3:end)), sub_MF_pt);
taskPerformanceDataPath = fullfile(taskperf_datadir, 'TaskPerformanceSD.mat');
load(taskPerformanceDataPath);
[~,~,IB_taskperformance] = intersect(sub_MF_pt_num, accuracyResults.Subject);
accuracyResults_updt = accuracyResults(IB_taskperformance,:);
meanAccBySubject_Animal = groupsummary(accuracyResults_updt, 'Subject', 'mean', 'Animal_ACC');
meanAccBySubject_Falsefont = groupsummary(accuracyResults_updt, 'Subject', 'mean', 'Falsefont_ACC');
totalmean = mean(meanAccBySubject_Falsefont.mean_Falsefont_ACC);
disp(['The total mean of mean_Falsefont_ACC is: ', num2str(totalmean)]);
totalstd = std(meanAccBySubject_Falsefont.mean_Falsefont_ACC);
disp(['The total std of mean_Falsefont_ACC is: ', num2str(totalstd)]);
totalmean = mean(meanAccBySubject_Animal.mean_Animal_ACC);
disp(['The total mean of mean_Anim_ACC is: ', num2str(totalmean)]);
totalstd = std(meanAccBySubject_Animal.mean_Animal_ACC);
disp(['The total std of mean_Falsefont_ACC is: ', num2str(totalstd)]);

%% Task Performance (2)
figure;
daboxplot(meanAccBySubject_Animal.mean_Animal_ACC, 'groups', ones(1, numel(meanAccBySubject_Animal.mean_Animal_ACC)), 'outsymbol', 'kx', 'xtlabels', 'test', 'fill', 0);
xlabel('All Subjects');
ylabel('Accuracy (%)');
set(gca, 'FontSize', 10);
hold on;
hold off;
title({'Task performance (anim)', 'of Discordant Subjects'});
set(gca, 'XTick', []);
box off;
set(gcf, 'Position', [1000   100   300   300]);
set(gca,'color','none');

figure;
daboxplot(meanAccBySubject_Falsefont.mean_Falsefont_ACC, 'groups', ones(1, numel(meanAccBySubject_Falsefont.mean_Falsefont_ACC)), 'outsymbol', 'kx', 'xtlabels', 'test', 'fill', 0);
xlabel('All Subjects');
ylabel('Accuracy (%)');
set(gca, 'FontSize', 10);
hold on;
hold off;
title({'Task performance (symb)', 'of Discordant Subjects'});
box off;
set(gcf, 'Position', [1000   100   300   300]);
set(gca,'color','none');

%% Corr (Response ACC, Reaction Time)
% Convert numerical subject IDs in accuracyResults_updt to string format with 'EC' prefix
accuracyResults_updt.SubjectStr = strcat('EC', arrayfun(@num2str, accuracyResults_updt.Subject, 'UniformOutput', false));

% Find the intersection between subject IDs
[commonSubIDs, idxAccuracy, idxRT] = intersect(accuracyResults_updt.SubjectStr, T_patn_MEGfMRI.Sub_ID);

% Extract relevant data
animalAcc = accuracyResults_updt.Animal_ACC(idxAccuracy);
falsefontAcc = accuracyResults_updt.Falsefont_ACC(idxAccuracy);
animalRT = T_patn_MEGfMRI.Animal(idxRT);
falsefontRT = T_patn_MEGfMRI.Symbol(idxRT);

% Calculate correlation for animal condition
[correlationAnimal, pValueAnimal] = corr(animalAcc, animalRT, 'Rows', 'complete');

% Calculate correlation for falsefont condition
[correlationFalsefont, pValueFalsefont] = corr(falsefontAcc, falsefontRT, 'Rows', 'complete');

% Display results
disp(['Correlation between Animal Accuracy and Reaction Time: ', num2str(correlationAnimal)]);
disp(['P-value for Animal Accuracy and Reaction Time: ', num2str(pValueAnimal)]);

disp(['Correlation between Falsefont Accuracy and Reaction Time: ', num2str(correlationFalsefont)]);
disp(['P-value for Falsefont Accuracy and Reaction Time: ', num2str(pValueFalsefont)]);

% Plot Correlations

% Plot for Animal condition
figure;
subplot 211
scatter(animalAcc, animalRT, 'filled');
title('Corr (Animal Acc, ReactionT');
xlabel('Animal Acc');
ylabel('ReactionT (s)');
grid on;
set(gca,'color','none');

% Plot for Falsefont condition
subplot 212
scatter(falsefontAcc, falsefontRT, 'filled');
title('Corr (Falsefont Acc, ReactionT)');
xlabel('Falsefont Acc');
ylabel('ReactionT (s)');
grid on;
box off;
set(gcf, 'Position', [1000   100   300   300]);
set(gca,'color','none');
