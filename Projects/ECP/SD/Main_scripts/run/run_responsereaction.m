% load('/data/MEG/Research/aizadi/process/RT_summary/ResponseTime.mat')
% [~,~,IB_reactiontime] = intersect(sub_MF_pt, T.Sub_ID);
% T_patn_MEGfMRI = T(IB_reactiontime,:);
% meanAnimal = mean(T_patn_MEGfMRI.Animal, 'omitnan');
% stdAnimal = std(T_patn_MEGfMRI.Animal, 'omitnan');
% meanSymbol = mean(T_patn_MEGfMRI.Symbol, 'omitnan');
% stdSymbol = std(T_patn_MEGfMRI.Symbol, 'omitnan');
% disp(['Mean of Animal reaction times: ', num2str(meanAnimal)]);
% disp(['Standard Deviation of Animal reaction times: ', num2str(stdAnimal)]);
% disp(['Mean of Symbol reaction times: ', num2str(meanSymbol)]);
% disp(['Standard Deviation of Symbol reaction times: ', num2str(stdSymbol)]);
% [correlationCoefficient, p] = corr(T_patn_MEGfMRI.Animal, T_patn_MEGfMRI.Symbol, 'Rows', 'complete');
% validPairs = sum(~isnan(T_patn_MEGfMRI.Animal) & ~isnan(T_patn_MEGfMRI.Symbol));
% df = validPairs - 2;
% disp(['Correlation coefficient: ', num2str(correlationCoefficient)]);
% disp(['Degrees of freedom: ', num2str(df)]);
% T_patn_MEGfMRI.Properties.VariableNames{'Sub_ID'} = 'SubjectID';
% T_patn_MEGfMRI.Properties.VariableNames{'Animal'} = 'Animal_RT';
% T_patn_MEGfMRI.Properties.VariableNames{'Symbol'} = 'Symbol_RT';




% load('/data/MEG/Research/ECP/Semantic_Decision/process/RT_summary/ResponseTime_meanRuns.mat'); 
% [~,~,IB_reactiontime] = intersect(sub_MF_pt, rt.sub_ID);
% length(IB_reactiontime)
% 
% T_patn_MEGfMRI = [];
% T_patn_MEGfMRI.Sub_ID = rt.sub_ID(IB_reactiontime);
% T_patn_MEGfMRI.Animal = rt.animal(IB_reactiontime);
% T_patn_MEGfMRI.Symbol = rt.symbol(IB_reactiontime);
% 
% % T_patn_MEGfMRI = T(IB_reactiontime,:);
% meanAnimal = mean(T_patn_MEGfMRI.Animal, 'omitnan');
% stdAnimal = std(T_patn_MEGfMRI.Animal, 'omitnan');
% meanSymbol = mean(T_patn_MEGfMRI.Symbol, 'omitnan');
% stdSymbol = std(T_patn_MEGfMRI.Symbol, 'omitnan');
% disp(['Mean of Animal reaction times: ', num2str(meanAnimal)]);
% disp(['Standard Deviation of Animal reaction times: ', num2str(stdAnimal)]);
% disp(['Mean of Symbol reaction times: ', num2str(meanSymbol)]);
% disp(['Standard Deviation of Symbol reaction times: ', num2str(stdSymbol)]);
% [correlationCoefficient, p] = corr(T_patn_MEGfMRI.Animal', T_patn_MEGfMRI.Symbol', 'Rows', 'complete');
% validPairs = sum(~isnan(T_patn_MEGfMRI.Animal) & ~isnan(T_patn_MEGfMRI.Symbol));
% df = validPairs - 2;
% disp(['Correlation coefficient: ', num2str(correlationCoefficient)]);
% disp(['Degrees of freedom: ', num2str(df)]);
% 
% T_patn_MEGfMRI = struct2table(T_patn_MEGfMRI);
% 
% T_patn_MEGfMRI.Properties.VariableNames{'Sub_ID'} = 'SubjectID';
% T_patn_MEGfMRI.Properties.VariableNames{'Animal'} = 'Animal_RT';
% T_patn_MEGfMRI.Properties.VariableNames{'Symbol'} = 'Symbol_RT';
% 
% disp(T_patn_MEGfMRI)

%%
% Load run-averaged RTs (assumes struct 'rt' with fields: sub_ID, animal, symbol, [both])
load('/data/MEG/Research/ECP/Semantic_Decision/process/RT_summary/ResponseTime_meanRuns.mat');

% Ensure compatible types
sub_MF_pt = string(sub_MF_pt(:));
rt_sub     = string(rt.SubjectID(:));
rt_animal  = rt.animal(:);
rt_symbol  = rt.symbol(:);
rt_both    = [];
hasBoth    = isfield(rt,'both');
if hasBoth, rt_both = rt.both(:); end

% Align to sub_MF_pt order and keep only matches
[commonSubs, ia, ib] = intersect(sub_MF_pt, rt_sub, 'stable');
fprintf('Matched subjects: %d\n', numel(ib));

% Build a proper row-per-subject table
Animal_RT = rt_animal(ib);
Symbol_RT = rt_symbol(ib);
T_patn_MEGfMRI = table(commonSubs, Animal_RT, Symbol_RT, 'VariableNames', ...
    {'SubjectID','Animal_RT','Symbol_RT'});

% Optionally include 'both' if present
if hasBoth
    Both_RT = rt_both(ib);
    T_patn_MEGfMRI.Both_RT = Both_RT;
end

% Summary stats
meanAnimal = mean(T_patn_MEGfMRI.Animal_RT,'omitnan');
stdAnimal  = std( T_patn_MEGfMRI.Animal_RT,'omitnan');
meanSymbol = mean(T_patn_MEGfMRI.Symbol_RT,'omitnan');
stdSymbol  = std( T_patn_MEGfMRI.Symbol_RT,'omitnan');
disp(['Mean Animal RT: ', num2str(meanAnimal), '  SD: ', num2str(stdAnimal)]);
disp(['Mean Symbol RT: ', num2str(meanSymbol), '  SD: ', num2str(stdSymbol)]);

% Correlation (complete cases)
[correlationCoefficient, p] = corr(T_patn_MEGfMRI.Animal_RT, T_patn_MEGfMRI.Symbol_RT, 'Rows','complete');
validPairs = sum(all(~isnan(T_patn_MEGfMRI{:,{'Animal_RT','Symbol_RT'}}),2));
df = validPairs - 2;
disp(['Correlation coefficient: ', num2str(correlationCoefficient), ', p = ', num2str(p)]);
disp(['Degrees of freedom: ', num2str(df)]);

% Inspect table
disp(T_patn_MEGfMRI)

% (Optional) save
% writetable(T_patn_MEGfMRI, '/data/MEG/Research/ECP/Semantic_Decision/process/RT_summary/RT_MEGfMRI_table.csv');
