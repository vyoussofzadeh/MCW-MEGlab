function data_in = SM_adjust_time(data_in,epoch_time)

% l = length(data_in.trial{1});
data_in.time = [];
for i = 1:length(data_in.trial)
    data_in.time{i} = linspace(0,epoch_time,epoch_time*data_in.fsample);
    data_in.trial{i} = data_in.trial{i}(:,1:epoch_time*data_in.fsample);
end