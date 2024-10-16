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