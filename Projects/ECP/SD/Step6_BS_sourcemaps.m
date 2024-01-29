

% Assuming your data is stored in a variable named 'brainstormData'
startTime = 0; % Example start time in seconds
endTime = 0.2;   % Example end time in seconds

[avgTMap, avgPMap] = calculateIntervalAverage(brainstormData, startTime, endTime);

% tmp = avgTMap;
% tmp = (tmp - min(tmp(:))) ./ (max(tmp(:)) - min(tmp(:)));
% avgTMap = tmp;

% avgTMap (avgTMap < 0.7.*max(avgTMap(:))) = 0;

newbrainstormData = brainstormData;

newbrainstormData.tmap = avgTMap;  % Replace 'data' with the actual variable name if different
newbrainstormData.pmap = avgPMap;

newbrainstormData.Comment = ['avg_', num2str(startTime), '_', num2str(endTime)]

newbrainstormData.Time = [1 1];