
function newbrainstormData = calculateIntervalWeightedAverage(brainstormData, startTime, endTime)
    % Find indices corresponding to the time interval
    timeIndices = brainstormData.Time >= startTime & brainstormData.Time <= endTime;

    % Extract data interval for source averages
    avgSourceData = mean(brainstormData.ImageGridAmp(:, timeIndices), 2);

    % Assuming you want to calculate avgTMap and avgPMap similarly
    % Update this part if the calculation for pmap is different
    avgTMap = avgSourceData;
    avgPMap = avgTMap;

    % Update brainstormData struct
    newbrainstormData = brainstormData;
    newbrainstormData.ImageGridAmp = avgSourceData; % Averaged source data
    newbrainstormData.tmap = avgTMap;  % Average tmap data
    newbrainstormData.pmap = avgPMap;  % Average pmap data
    newbrainstormData.Comment = ['Interval Avg: ', num2str(startTime), ' to ', num2str(endTime)];
    newbrainstormData.Time = [startTime endTime]; % Updated time to reflect the averaging interval
end
