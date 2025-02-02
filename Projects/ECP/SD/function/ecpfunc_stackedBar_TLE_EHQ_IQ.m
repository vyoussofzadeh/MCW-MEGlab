function ecpfunc_stackedBar_TLE_EHQ_IQ(bestResultsTable, T1_epil_measures)
% ECPFUNC_STACKEDBAR_TLE_EHQ_IQ
%
% Creates:
%   1) A stacked bar of TLE side (Left, Bilateral, Right) 
%      for discordant subjects in each ROI.
%   2) Converts numeric EHQ (-100..+100) to categorical 
%      (Left/Ambi/Right) using ±40 thresholds, then shows 
%      a stacked bar by handedness for discordant subjects in each ROI.
%   3) A box plot of IQ (e.g., NP1WASI_FSIQ) for discordant subjects 
%      in each ROI, as an example of how to add more continuous measures.
%
% USAGE:
%   ecpfunc_stackedBar_TLE_EHQ_IQ(bestResultsTable, T1_epil_measures);
%
% INPUTS:
%   bestResultsTable : A table with at least the columns:
%       - ROI              : ROI names (cell array of strings)
%       - Best_Discord_Subs: cell array of numeric indices of discordant subs
%   T1_epil_measures : A table with at least:
%       - TLEside  : Nx1 categorical ('Left','Right','Bilateral')
%       - EHQ      : Nx1 numeric (-100..+100)
%       - NP1WASI_FSIQ : Nx1 numeric (IQ measure)
%
% Author: ChatGPT, adapted from your code snippets (2025).

%% 1) Basic check
if ~istable(bestResultsTable) || ~istable(T1_epil_measures)
    error('Inputs must be tables.');
end
if ~all(ismember({'ROI','Best_Discord_Subs'}, bestResultsTable.Properties.VariableNames))
    error('bestResultsTable must have columns ROI, Best_Discord_Subs.');
end
if ~all(ismember({'TLEside','EHQ','NP1WASI_FSIQ'}, T1_epil_measures.Properties.VariableNames))
    error('T1_epil_measures must have columns TLEside, EHQ, NP1WASI_FSIQ.');
end

%% 2) TLE side stacked bar (Left / Bilateral / Right) among discordant
allROI = bestResultsTable.ROI;     % cell array of ROI names
nROIs  = height(bestResultsTable); % number of ROIs

% Pre-allocate counters
nLeftTLE  = zeros(nROIs,1);
nRightTLE = zeros(nROIs,1);
nBilatTLE = zeros(nROIs,1);

for i = 1:nROIs
    
    % A) discordant subject indices for this ROI
    discordSubs = bestResultsTable.Best_Discord_Subs{i};
    
    % B) TLE side for these subjects
    TLEside_discord = T1_epil_measures.TLEside(discordSubs);
    
    % C) Count how many have Left vs Bilateral vs Right
    nLeftTLE(i)   = sum(TLEside_discord == 'Left');
    nBilatTLE(i)  = sum(TLEside_discord == 'Bilateral');
    nRightTLE(i)  = sum(TLEside_discord == 'Right');
end

% D) Plot a stacked bar
figure('Color','w','Name','Discordant By TLE Side','Position',[100 300 700 400]);
barData = [nLeftTLE, nBilatTLE, nRightTLE];
hBar = bar(barData, 'stacked');

% Optional: color each stack
hBar(1).FaceColor = [0.8 0.2 0.2];   % red
hBar(2).FaceColor = [0.5 0.5 0.5];   % gray
hBar(3).FaceColor = [0.2 0.6 0.8];   % teal

set(gca, 'XTick', 1:nROIs, 'XTickLabel', allROI, 'FontSize', 10);
xlabel('ROI');
ylabel('# Discordant Patients');
title('Discordant Subjects by TLE Side');
legend({'Left TLE','Bilateral','Right TLE'}, 'Location','bestoutside');
box off;

%% 3) Convert numeric EHQ to categorical (Left/Ambi/Right) using ±40 thresholds
ehqVals  = T1_epil_measures.EHQ;
ehqCatStr = repmat("Ambi", size(ehqVals));  % default
threshold = 40;

ehqCatStr(ehqVals >  threshold) = "Right";
ehqCatStr(ehqVals < -threshold) = "Left";

T1_epil_measures.EHQcat = categorical(ehqCatStr, ["Left","Ambi","Right"]);

%% 4) EHQ stacked bar (Left-Handed, Ambi, Right-Handed) among discordant
nLeftEHQ  = zeros(nROIs,1);
nAmbiEHQ  = zeros(nROIs,1);
nRightEHQ = zeros(nROIs,1);

for i = 1:nROIs
    discordSubs = bestResultsTable.Best_Discord_Subs{i};
    ehqCat_discord = T1_epil_measures.EHQcat(discordSubs);

    nLeftEHQ(i)  = sum(ehqCat_discord == 'Left');
    nAmbiEHQ(i)  = sum(ehqCat_discord == 'Ambi');
    nRightEHQ(i) = sum(ehqCat_discord == 'Right');
end

figure('Color','w','Name','Discordant By EHQ Category','Position',[100 300 800 400]);
barDataEHQ = [nLeftEHQ, nAmbiEHQ, nRightEHQ];
hBar2 = bar(barDataEHQ, 'stacked');

% Optional colors
hBar2(1).FaceColor = [0.8 0.2 0.2];   % red
hBar2(2).FaceColor = [0.5 0.5 0.5];   % gray
hBar2(3).FaceColor = [0.2 0.6 0.8];   % teal

set(gca, 'XTick', 1:nROIs, 'XTickLabel', allROI, 'FontSize', 10);
xlabel('ROI');
ylabel('# Discordant Subjects');
title('Discordant Subjects by Handedness (EHQ)');
legend({'Left-Handed','Ambidextrous','Right-Handed'}, 'Location','bestoutside');
box off;

%% 5) Example: Box plot of IQ (NP1WASI_FSIQ) for discordant subs across ROIs
%    We'll gather each ROI's discordant subs' IQ, then do a grouped box plot.

allIQ = {};     % cell array, each element is a vector of IQs for that ROI
for i = 1:nROIs
    discordSubs = bestResultsTable.Best_Discord_Subs{i};
    roiIQ = T1_epil_measures.NP1WASI_FSIQ(discordSubs);
    allIQ{i} = roiIQ;
end

% We can combine them into one numeric array and provide a grouping vector
allIQvalues = []; 
group = [];       % which ROI each data point belongs to
for i = 1:nROIs
    allIQvalues = [allIQvalues; allIQ{i}];
    group = [group; repmat(i, length(allIQ{i}), 1)];
end

% Now boxplot with 'group' defining each ROI
figure('Color','w','Name','IQ Box Plot','Position',[100 300 600 400]);
boxplot(allIQvalues, group, 'Labels', allROI, 'Whisker',1.5);

xlabel('ROI');
ylabel('IQ (NP1WASI-FSIQ)');
title('Box Plot of IQ for Discordant Subjects by ROI');
box off;

end
