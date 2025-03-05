%% MEG noise-SNR
% Load Data and Compute med_snr_mag / med_snr_grad
% cd('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/data/SNR_Data')
cd('/data/MEG/Research/MEGIN project/Scripts/SNR/SNR_Data')
[all_snr_grad, all_snr_mag] = ecpfunc_read_noiseSNR();

mean_snr_mag  = (all_snr_mag(1).snr_value  + all_snr_mag(2).snr_value)  ./ 2;  
mean_snr_grad = (all_snr_grad(1).snr_value + all_snr_grad(2).snr_value) ./ 2;

mean_snr_MEGnetvstSSS  = [all_snr_mag(1).snr_value,  all_snr_grad(1).snr_value];  
mean_snr_tSSSvsRaw = [all_snr_mag(2).snr_value, all_snr_grad(2).snr_value];

size(mean_snr_tSSSvsRaw)
size(mean_snr_MEGnetvstSSS)

% med_snr_mag and med_snr_grad => [N_subj x N_chan]
med_snr_mag  = squeeze(median(mean_snr_mag, 2));   
med_snr_grad = squeeze(median(mean_snr_grad, 2)); 

med_snr_tSSSvsRaw  = squeeze(median(mean_snr_MEGnetvstSSS, 2));   
med_snr_MEGnetvstSSS = squeeze(median(mean_snr_tSSSvsRaw, 2)); 

% "mean_mag" => average across channels => [N_subj x 1]
mean_mag  = mean(med_snr_mag,  2); 
mean_grad = mean(med_snr_grad, 2);

mean_tSSSvsRaw  = mean(med_snr_tSSSvsRaw,  2); 
mean_MEGnetvstSSS = mean(med_snr_MEGnetvstSSS, 2);

% If you also want a combined metric for each subject (e.g. sum or average):
% Here, for example, we add them and then take the mean across channels:
%   mean_snr = mean(med_snr_mag + med_snr_grad, 2);
% But if you've already used "mean(...,2)" above, just do:
mean_snr = (mean_mag + mean_grad) / 2;
mean_snr2 = (mean_tSSSvsRaw + mean_MEGnetvstSSS) / 2;

nSubjects = length(mean_snr);

% Identify or retrieve subject IDs:
sub_snr_list = unique(all_snr_mag(1).sub);  % cell array of length N_subj
if length(sub_snr_list) ~= nSubjects
    error('Mismatch: sub_snr_list has %d entries, but mean_snr has %d.', ...
          length(sub_snr_list), nSubjects);
end

figure;
% bar(mean_grad);
scatter(1:nSubjects, mean_grad, 'filled');
title('Mean SNR Across Subjects');
xlabel('Subjects');
ylabel('SNR (absolute value)');
xticks(1:nSubjects);
xticklabels(sub_snr_list);
xtickangle(90);
set(gca, 'FontSize', 8);

figure;
% bar(mean_mag);
scatter(1:nSubjects, mean_mag, 'filled');
title('Mean SNR Across Subjects');
xlabel('Subjects');
ylabel('SNR (absolute value)');
xticks(1:nSubjects);
xticklabels(sub_snr_list);
xtickangle(90);
set(gca, 'FontSize', 8);

figure;
% bar(mean_grad);
scatter(1:nSubjects, mean_tSSSvsRaw, 'filled');
title('Mean SNR Across Subjects, mean_tSSSvsRaw');
xlabel('Subjects');
ylabel('SNR (absolute value)');
xticks(1:nSubjects);
xticklabels(sub_snr_list);
xtickangle(90);
set(gca, 'FontSize', 8);

figure;
% bar(mean_mag);
scatter(1:nSubjects, mean_MEGnetvstSSS, 'filled');
title('Mean SNR Across Subjects, mean_MEGnetvstSSS');
xlabel('Subjects');
ylabel('SNR (absolute value)');
xticks(1:nSubjects);
xticklabels(sub_snr_list);
xtickangle(90);
set(gca, 'FontSize', 8);