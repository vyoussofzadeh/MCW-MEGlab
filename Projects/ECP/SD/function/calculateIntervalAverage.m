% function [avgTMap, avgPMap] = calculateIntervalAverage(data, startTime, endTime)
%     % Find indices corresponding to the time interval
%     timeIndices = find(data.Time >= startTime & data.Time <= endTime);
%
%     % Extract the relevant parts of the tmap and pmap
%     tMapInterval = data.tmap(:, timeIndices);
%     pMapInterval = data.pmap(:, timeIndices);
%
%     % Calculate the average across the time interval
%     avgTMap = mean(tMapInterval, 2);
%     avgPMap = mean(pMapInterval, 2);
% end


function newbrainstormData = calculateIntervalAverage(brainstormData, startTime, endTime)
    % Find indices corresponding to the time interval
    timeIndices = brainstormData.Time >= startTime & brainstormData.Time <= endTime;

    % Extract data interval for source averages
%     avgSourceData = mean(brainstormData.ImageGridAmp(:, timeIndices), 2);
    avgpmap = mean(brainstormData.pmap(:, timeIndices), 2);
    avgtmap = mean(brainstormData.tmap(:, timeIndices), 2);


    % Assuming you want to calculate avgTMap and avgPMap similarly
    % Update this part if the calculation for pmap is different
%     avgTMap = avgSourceData;
    avgTMap = avgtmap;
    avgPMap = avgpmap;

    % Update brainstormData struct
    newbrainstormData = brainstormData;
%     newbrainstormData.ImageGridAmp = avgTMap; % Averaged source data
    newbrainstormData.tmap = avgTMap;  % Average tmap data
    newbrainstormData.pmap = avgPMap;  % Average pmap data
    newbrainstormData.Comment = ['Interval Avg: ', num2str(startTime), ' to ', num2str(endTime)];
    newbrainstormData.Time = [startTime endTime]; % Updated time to reflect the averaging interval
end


% function newbrainstormData = calculateIntervalAverage(brainstormData, startTime, endTime)
% % Find indices corresponding to the time interval
% timeIndices = brainstormData.Time >= startTime & brainstormData.Time <= endTime;
% 
% % Assuming ImageGridAmp contains the data of interest
% % You might need to adjust this if your actual tmap and pmap data are located elsewhere
% dataInterval = brainstormData.ImageGridAmp(:, timeIndices);
% 
% % Calculate the average across the time interval
% % Modify this part if you need to calculate avgTMap and avgPMap separately
% % Here, I'm assuming you're averaging the same data, but the structure of your data might require different handling
% avgTMap = mean(dataInterval, 2);
% avgPMap = avgTMap; % If pmap calculation is different, update this line accordingly
% 
% newbrainstormData = brainstormData;
% 
% newbrainstormData.tmap = avgTMap;  % Replace 'data' with the actual variable name if different
% newbrainstormData.pmap = avgPMap;
% 
% newbrainstormData.Comment = ['avg_', num2str(startTime), '_', num2str(endTime)];
% 
% newbrainstormData.Time = [1 1];
% end

