function dataCell = do_normalize_data(dataCell)


for i = 1:size(dataCell,1)
    dataCell(i,:) = (dataCell(i,:) - min(dataCell(i,:), [], 'all')) / (max(dataCell(i,:), [], 'all') - min(dataCell(i,:), [], 'all'));
end

