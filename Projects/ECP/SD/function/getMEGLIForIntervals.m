function MEG_LI_cell = getMEGLIForIntervals(MEG_LI, summaryTable, wi)
% GETMEGLIFORINTERVALS Extracts MEG_LI columns for each time interval 
% specified in summaryTable, based on wi.
%
%   MEG_LI_cell = getMEGLIForIntervals(MEG_LI, summaryTable, wi)
%
%   INPUTS:
%       MEG_LI       - (#Subjects x #Intervals) numeric matrix of MEG LI values
%                      where columns correspond to rows in wi (e.g., 44 intervals).
%       summaryTable - A table with a column named 'Time_Interval'
%                      containing [start end], e.g. 0.3  0.6
%       wi           - (#Intervals x 2) numeric matrix of interval boundaries,
%                      where wi(k,:) = [start_k, end_k] is the k-th interval.
%
%   OUTPUT:
%       MEG_LI_cell  - cell array of size (height(summaryTable) x 1),
%                      each cell containing the sub-matrix of MEG_LI 
%                      (rows=#Subjects, columns=the matched intervals).
%
%   EXAMPLE:
%       % Suppose MEG_LI is 72 x 44, summaryTable has a 'Time_Interval' column, 
%       % and wi is 44 x 2. Then:
%       MEG_LI_cell = getMEGLIForIntervals(MEG_LI, summaryTable, wi);
%       % MEG_LI_cell{i} is the (#Subjects x #Cols) submatrix for row i in summaryTable.
%
%   Author: (Your Name / Organization)
%   Date: (Date)

    % Validate that summaryTable has a 'Time_Interval' column
    if ~ismember('Time_Interval', summaryTable.Properties.VariableNames)
        error('summaryTable must have a "Time_Interval" column with [start end].');
    end

    nRows = height(summaryTable);
    MEG_LI_cell = cell(nRows,1);

    % For each row of summaryTable
    for iRow = 1:nRows
        % Extract the [startTime endTime] from summaryTable.Time_Interval
        intervalVal = summaryTable.Time_Interval(iRow, :);
        startTime   = intervalVal(1);
        endTime     = intervalVal(2);

        % (A) EXACT MATCH scenario:
        % If you want intervals that exactly match the row in wi:
        % colMask = (wi(:,1) == startTime & wi(:,2) == endTime);
        %
        % (B) RANGE scenario:
        % If you want intervals whose [wi(k,1) wi(k,2)] is fully inside [startTime endTime]:
        % colMask = (wi(:,1) >= startTime & wi(:,2) <= endTime);
        %
        % (C) Or you can interpret partial overlap, etc.
        % For a simple approach, let's assume EXACT MATCH:

        colMask = (wi(:,1) == startTime & wi(:,2) == endTime);

        % Now, find which columns satisfy colMask
        matchedCols = find(colMask);
%         matchedCols_all(iRow) = matchedCols;

        if isempty(matchedCols)
            % No exact match found in wi for [startTime, endTime]
            % We'll leave MEG_LI_cell{iRow} = [] or store NaN
            MEG_LI_cell{iRow} = [];
        else
            % Extract those columns from MEG_LI
            MEG_LI_cell{iRow} = MEG_LI(:, matchedCols);
        end
    end
end
