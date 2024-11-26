function triggerChannels = detectTriggerChannels(data)

evt = ft_read_event(data);
hdr = ft_read_header(data); %read header information

types = arrayfun(@(x) x.type, evt, 'UniformOutput', false);
% Finding unique 'type' values
unique_types = unique(types);

% Assume 'data' is a struct with fields 'times' and 'channels'
triggerChannels = struct(); % Initialize output struct

for i = 1:length(unique_types)
    channelName = unique_types{i};
    if startsWith(channelName, 'STI') % Check only STI* channels
        channelData = data.channels.(channelName);
        % A simple trigger detection strategy:
        % Check if there are any non-zero values and variability
        if any(channelData ~= 0) && numel(unique(channelData)) > 1
            fprintf('Trigger data detected in %s\n', channelName);
            triggerChannels.(channelName) = channelData;
        end
    end
end
end
