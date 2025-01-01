function do_combine_images(folderPath, outputFile)
% do_combine_images Combines images from a specified folder into one single image 
% with auto-resizing, a white background, and optional manual re-ordering.
%
% Usage:
%   do_combine_images(folderPath, outputFile)
%
% Inputs:
%   folderPath - (String) Path to the folder containing images (e.g., PNGs).
%   outputFile - (String) Name of the final merged file (e.g., 'combined.png').
%
% The function lists all images found, gives you an opportunity to
% reorder them (by typing a custom sequence of indices), and merges
% them into a grid saved in outputFile.

    % 1. Read image files
    files = dir(fullfile(folderPath, '*.png'));  % Adjust extension if needed
    if isempty(files)
        error('No images found in the specified folder.');
    end

    % 2. Sort files by name (alphabetically), then exclude outputFile
    [~, idx] = sort({files.name});
    files = files(idx);
    files = files(~strcmp({files.name}, outputFile));

    % 3. Sort remaining files by modification date (ascending)
    fileDates = [files.datenum];
    [~, idx] = sort(fileDates);
    files = files(idx);
    
    % 4. Display the files to be merged
    fprintf('The following images will be merged (current order):\n');
    for i = 1:length(files)
        fprintf('%2d: %s\n', i, files(i).name);
    end
    
    % 5. Ask user for a custom re-ordering
%     promptMsg = sprintf( ...
%         '\nEnter a comma-separated list of indices (e.g. 3,1,2,...) to re-order,\n' ...
%         'or press Enter to keep the current order: ');
%     userInput = input(promptMsg, 's');
    
    promptMsg = sprintf('\nEnter a comma-separated list of indices (e.g. 3,1,2,...) to re-order,\nor press Enter to keep the current order: ');
    
    userInput = input(promptMsg, 's');
    
    
    if ~isempty(userInput)
        % Parse user input into numeric array
        newOrder = str2num(userInput); %#ok<ST2NM> 
        if isempty(newOrder) || length(newOrder) ~= length(files) ...
                || any(newOrder < 1) || any(newOrder > length(files)) ...
                || length(unique(newOrder)) ~= length(files)
            error(['Invalid indices. Ensure you provide a comma-separated list ' ...
                   'with each index from 1 to %d exactly once.'], length(files));
        end
        files = files(newOrder);
        fprintf('Images re-ordered by user input.\n');
    else
        fprintf('Proceeding with default order.\n');
    end

    % 6. Print final order
    fprintf('\nFinal merge order:\n');
    for i = 1:length(files)
        fprintf('%2d: %s\n', i, files(i).name);
    end
    fprintf('\n');

    % 7. Merge images into a grid
    numImages = length(files);
    numPerRow = ceil(sqrt(numImages));  % # of images per row/column

    % Read first image to get reference dimensions
    firstImage = imread(fullfile(folderPath, files(1).name));
    [rows, cols, channels] = size(firstImage);

    % Initialize white background
    if isa(firstImage, 'uint8')
        combinedImage = 255 * ones(rows*numPerRow, cols*numPerRow, channels, 'uint8');
    else
        combinedImage = 255 * ones(rows*numPerRow, cols*numPerRow, channels, 'like', firstImage);
    end

    % 8. Place each image in the combined grid
    for k = 1:numImages
        % Grid position (row/col in the final big image)
        gridRow = floor((k-1) / numPerRow);
        gridCol = mod(k-1, numPerRow);

        % Read and resize the image to match reference dimensions
        thisImage = imread(fullfile(folderPath, files(k).name));
        resizedImage = imresize(thisImage, [rows, cols]);

        % Insert it into combinedImage
        rowStart = gridRow*rows + 1;
        colStart = gridCol*cols + 1;
        combinedImage(rowStart:(rowStart+rows-1), ...
                      colStart:(colStart+cols-1), :) = resizedImage;
    end

    % 9. Write the combined image
    export_path = fullfile(folderPath, outputFile);
    imwrite(combinedImage, export_path);

    % 10. Confirm
    fprintf('Combined image saved as: %s\n', export_path);
    cd(folderPath);
end


% function do_combine_images(folderPath, outputFile)
% % Combine images from a specified folder into one single image with auto-resizing and a white background.
% 
% % Parameters:
% % folderPath - String, path to the folder containing images.
% % outputFile - String, path and filename of the output image.
% 
% % Read all image files from the folder
% files = dir(fullfile(folderPath, '*.png'));  % Adjust the extension if needed
% if isempty(files)
%     error('No images found in the specified folder.');
% end
% 
% % Sort files by name to ensure they are processed in alphabetical order
% [~, idx] = sort({files.name});
% files = files(idx);
% 
% % Filter out the outputFile if it exists in the same folder
% files = files(~strcmp({files.name}, outputFile));
% 
% % Extract date info and sort files by date
% fileDates = [files.datenum];
% [~, idx] = sort(fileDates);
% files = files(idx);
% 
% % Assuming you've already sorted `files` as needed
% fprintf('The following images will be merged:\n');
% for i = 1:length(files)
%     fprintf('%d: %s\n', i, files(i).name);
% end
% 
% % Determine number of images
% numImages = length(files);
% numPerRow = ceil(sqrt(numImages));  % Number of images per row and column
% 
% % Read the first image to get dimensions and type
% firstImage = imread(fullfile(folderPath, files(1).name));
% [rows, cols, channels] = size(firstImage);
% 
% % Initialize the combined image with a white background
% if isa(firstImage, 'uint8')
%     combinedImage = 255 * ones(rows * numPerRow, cols * numPerRow, channels, 'uint8');
% else
%     combinedImage = 255 * ones(rows * numPerRow, cols * numPerRow, channels, 'like', firstImage);
% end
% 
% % Place each image into the combined image
% for k = 1:numImages
%     % Determine position in grid
%     row = floor((k-1) / numPerRow);
%     col = mod(k-1, numPerRow);
%     
%     % Read and resize image
%     thisImage = imread(fullfile(folderPath, files(k).name));
%     resizedImage = imresize(thisImage, [rows, cols]);  % Resize the image
%     
%     % Place image in the corresponding part of the combined image
%     combinedImage(row*rows + (1:rows), col*cols + (1:cols), :) = resizedImage;
% end
% 
% export_path = fullfile(folderPath,outputFile);
% 
% % Write the combined image to file
% imwrite(combinedImage, export_path);
% 
% % imshow(combinedImage)
% % print(gcf, export_path, '-dpng', '-r300')
% % print(gcf, export_path, '-depsc', '-r300')
% 
% fprintf('Combined image saved as %s\n', outputFile);
% cd(folderPath)
% end
