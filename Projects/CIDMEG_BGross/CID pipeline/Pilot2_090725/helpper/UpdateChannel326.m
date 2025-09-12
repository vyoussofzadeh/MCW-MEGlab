
%%
% Copy the structure from an existing channel (assuming index 1 is safe)
template = Chan324.Channel(1);

% Clear all fields (optional, to avoid confusion)
fields = fieldnames(template);
for f = 1:numel(fields)
    template.(fields{f}) = [];
end

% Create first dummy channel
dummy1 = template;
dummy1.Type = 'Misc';
dummy1.Name = 'MISC001';

% Create second dummy channel
dummy2 = template;
dummy2.Type = 'Misc';
dummy2.Name = 'MISC002';

% Append to existing channel list
Chan324.Channel(end+1) = dummy1;
Chan324.Channel(end+1) = dummy2;


Chan324.Comment = 'Neuromag channels (326)';