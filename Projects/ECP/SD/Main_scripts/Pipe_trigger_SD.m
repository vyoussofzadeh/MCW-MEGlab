clear;
clc;
close all;

ft_path = '/opt/matlab_toolboxes/ft_packages/fieldtrip_latest';
addpath(ft_path);
ft_defaults

% Define file path for the MEG data
sd_fif = '/group/jbinder/ECP/MEG/MEG_Work/EC1002/tSSS/tSSS_ec1002_SD_run1_raw.fif';

% Configuration for defining trials
cfg = [];
cfg.dataset = sd_fif;
cfg.trialfun = 'ft_trialfun_general'; % Default trial function
cfg.trialdef.eventtype = 'STI101'; % Event type for the stimulus trigger
% cfg.trialdef.eventvalue = cfg_main.eventid; % Uncomment and specify if needed
cfg.trialdef.prestim = 1; % Pre-stimulus time in seconds
cfg.trialdef.poststim = 3; % Post-stimulus time in seconds

% Define trials based on the configuration
evt = ft_definetrial(cfg);

% Define pre and post stimulus times
prestimTime = 1; 
poststimTime = 3;

% Read header information from the MEG file
hdr = ft_read_header(sd_fif);
Fsample = hdr.Fs; % Sampling frequency
prestimSamples = floor(prestimTime * Fsample); % Pre-stimulus samples
poststimSamples = floor(poststimTime * Fsample); % Post-stimulus samples

% Find the index of the trigger channel
Index = find(strcmp(hdr.label, 'STI101'));

% Read trigger data from the MEG file
detTrig = ft_read_data(sd_fif, 'chanindx', Index, 'header', hdr, 'eventformat', 'neuromag_fif', 'dataformat', 'neuromag_fif');
detTrig = bitand(full(detTrig),255);
% figure,plot(detTrig)

detResp=ft_read_data(sd_fif,'chanindx',Index,'header',hdr,'eventformat','digital trigger','dataformat','neuromag_fif');
detResp = detResp - detTrig;
% figure,plot(detResp)
Nsamps = length(detTrig);
% event = ft_read_event(datafile);

figure,
subplot 211
% plot(detTrig)
% hold on,
plot(detResp,'r')
legend({'response'})
title('Trigger Data');
xlabel('Time (samples)');
ylabel('Trigger value');

% Plot the trigger data
% figure;
subplot 212
plot(detTrig);
legend({'Trigger'});
title('Trigger Data');
xlabel('Time (samples)');
ylabel('Trigger value');


% Calculate the differences between consecutive samples
triggerDifferences = diff(detTrig);  % Pad with zero to align with original vector size


%%
% Calculate differences to find changes
triggerOnsets = find(triggerDifferences ~= 0); % Indices of changes
triggerTimes = triggerOnsets / Fsample;  % Convert indices to time in seconds
timeIntervals = diff(triggerTimes);  % Time intervals between triggers

% Display intervals
disp('Time intervals between triggers:');
disp(timeIntervals);

% Plot the intervals
figure;
stem(triggerTimes(1:end-1), timeIntervals);  % Plot intervals against their starting times
title('Time Intervals Between Trigger Events');
xlabel('Time (seconds)');
ylabel('Interval (seconds)');
grid on;

%%
val_idx = (timeIntervals < 1.5) & (timeIntervals > 0.5);

figure;
stem(timeIntervals(val_idx));  % Plot intervals against their starting times
title('Time Intervals Between Trigger Events');
xlabel('Time (seconds)');
ylabel('Interval (seconds)');
grid on;


triggerTimes_resp = triggerTimes(val_idx);

expectedWaitingTime = 1.0; % Expected interval in seconds

% Calculate jitter as the deviation from the expected interval
% Assuming 'triggerTimes' is a vector of times when each stimulus starts
timeIntervals = diff(triggerTimes);  % Time intervals between consecutive stimuli

% Adjust intervals to reflect only the waiting time
waitingTimes = timeIntervals; % Subtract the 2-second stimulus duration from the total intervals

% Calculate jitter as the deviation from the expected waiting time
jitterValues = waitingTimes - expectedWaitingTime;

% Display jitter values
disp('Jitter values for each waiting period:');
disp(jitterValues);

% Basic statistics
meanJitter = mean(jitterValues);
stdJitter = std(jitterValues);
maxJitter = max(jitterValues);
minJitter = min(jitterValues);

% Display results
disp(['Mean Jitter: ', num2str(meanJitter), ' seconds']);
disp(['Standard Deviation of Jitter: ', num2str(stdJitter), ' seconds']);
disp(['Maximum Jitter: ', num2str(maxJitter), ' seconds']);
disp(['Minimum Jitter: ', num2str(minJitter), ' seconds']);

figure;
histogram(jitterValues, 'BinWidth', 0.05); % Adjust the number of bins or bin width as needed
title('Distribution of Jitter Values');
xlabel('Jitter (seconds)');
ylabel('Frequency');
grid on;


%%

% Find indices where detTrig transitions to 3
% This assumes that the onset value '3' appears rising from zero or another baseline value.
sample_3 = find(triggerDifferences == 3);
sample_2 = find(triggerDifferences == 2);

% Combine and sort all sample indices
all_samples = sort([sample_3, sample_2]);
all_times = all_samples / Fsample; % Convert to time in seconds

% Calculate ITI for all trials
ITI_all = diff(all_times);

% Display mean ITI
mean_ITI_all = mean(ITI_all);
disp(['Mean ITI for all trials: ', num2str(mean_ITI_all), ' seconds']);
disp(['Range ITI for all trials: ', num2str(min(ITI_all)), '  , ', num2str(max(ITI_all))]);


% Compute the total silence time (ITI sum)
total_ITI_time = sum(ITI_all);
disp(['Total silence time between trials: ', num2str(total_ITI_time), ' seconds']);

% Compute and display timing of individual trials
disp('Timing of individual trials (seconds):');


% Compute and display timing of individual trials
disp('Timing of individual events and their respective silence times (seconds):');
for i = 1:length(all_samples) - 1
    
    trial_start = all_times(i);
    trial_end = all_times(i + 1);
    silence_time = ITI_all(i);

    % Display the start and end times of each trial and the silence period between them
    disp(['Trial ', num2str(i), ': Start = ', num2str(trial_start), ' s, End = ', num2str(trial_end), ' s, Silence Time = ', num2str(silence_time), ' s']);
end

% Plotting the trigger data over time with event markers
timeVector = (1:length(detTrig)) / Fsample; % Convert sample numbers to time in seconds

figure;
subplot(211);
plot(timeVector, detTrig);
hold on;
scatter(all_times, detTrig(all_samples), 'filled');  % Mark all sample indices on the plot
title('Trigger Data with Event Markers');
xlabel('Time (seconds)');
ylabel('Trigger Value');
legend('Trigger Signal', 'Detected Events');
hold off;

% Add a subplot for the ITI distribution over time
subplot(212);
stem(all_times(1:end-1), ITI_all, 'filled');
title('Inter-Trial Intervals (ITI) Over Time');
xlabel('Time (seconds)');
ylabel('ITI (seconds)');
grid on;

% Save the plot
savefig('Trigger_and_ITI_Plot.fig');

figure;
subplot 211
plot(timeVector, detTrig);
hold on;
scatter(sample_3 / Fsample, detTrig(sample_3), 'r', 'filled'); % Convert sample indices to time
scatter(sample_2 / Fsample, detTrig(sample_2), 'g', 'filled'); % Convert sample indices to time
ylim([-0.2,3.2])
lgnd = legend({'Trigger','Animal (3)', 'Symbol (2)'}, 'Location', 'south');
title('Trigger Data with Event Markers');
xlabel('Time (seconds)');
ylabel('Trigger value');
hold off;
set(gca,'color','none');
% set(lgnd,'color','none');

subplot 212
% plot(detTrig)
% hold on,
plot(timeVector, detResp,'r');
set(gca,'color','none');
legend({'response'})
title('Response Data');
xlabel('Time (samples)');
ylabel('Response value');

% Save the plots
addpath('/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/FT_fucntions/functions_new')
save_dir = '/data/MEG/Vahab/Github/MCW_MEGlab/MCW_MEGlab_git/Projects/ECP/SD/results';

cfg = []; 
cfg.outdir = save_dir; 
cfg.filename = 'SD_stim'; 
cfg.type = 'fig'; 
do_export_fig(cfg);

% Plotting ITI
figure;
stem(all_times(1:end-1), ITI_all, 'Marker', 'none');
title('Inter-Trial Intervals (ITI) over Time');
xlabel('Time (seconds)');
ylabel('ITI (seconds)');
grid on;

cfg = []; 
cfg.outdir = save_dir; 
cfg.filename = 'SD_ITI'; 
cfg.type = 'fig'; 
do_export_fig(cfg);


