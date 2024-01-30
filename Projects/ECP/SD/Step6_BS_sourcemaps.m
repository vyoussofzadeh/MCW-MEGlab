

% Assuming your data is stored in a variable named 'brainstormData'
startTime = -0.3; % Example start time in seconds
endTime = 0;   % Example end time in seconds

startTime = 0; % Example start time in seconds
endTime = .200;   % Example end time in seconds

startTime = .200; % Example start time in seconds
endTime = .400;   % Example end time in seconds

startTime = .400; % Example start time in seconds
endTime = .600;   % Example end time in seconds

startTime = .600; % Example start time in seconds
endTime = .800;   % Example end time in seconds

startTime = .800; % Example start time in seconds
endTime = 1;   % Example end time in seconds

newbrainstormData = calculateIntervalAverage(brainstormData, startTime, endTime);