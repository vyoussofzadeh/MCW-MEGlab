function do_combine_images(folderPath, outputFile)
    % Combine images from a specified folder into one single image with auto-resizing and a white background.
    
    % Parameters:
    % folderPath - String, path to the folder containing images.
    % outputFile - String, path and filename of the output image.
    
    % Read all image files from the folder
    files = dir(fullfile(folderPath, '*.png'));  % Adjust the extension if needed
    if isempty(files)
        error('No images found in the specified folder.');
    end
    
    % Sort files by name to ensure they are processed in alphabetical order
    [~, idx] = sort({files.name});
    files = files(idx);
    
    % Determine number of images
    numImages = length(files);
    numPerRow = ceil(sqrt(numImages));  % Number of images per row and column
    
    % Read the first image to get dimensions and type
    firstImage = imread(fullfile(folderPath, files(1).name));
    [rows, cols, channels] = size(firstImage);
    
    % Initialize the combined image with a white background
    if isa(firstImage, 'uint8')
        combinedImage = 255 * ones(rows * numPerRow, cols * numPerRow, channels, 'uint8');
    else
        combinedImage = 255 * ones(rows * numPerRow, cols * numPerRow, channels, 'like', firstImage);
    end
    
    % Place each image into the combined image
    for k = 1:numImages
        % Determine position in grid
        row = floor((k-1) / numPerRow);
        col = mod(k-1, numPerRow);
        
        % Read and resize image
        thisImage = imread(fullfile(folderPath, files(k).name));
        resizedImage = imresize(thisImage, [rows, cols]);  % Resize the image
        
        % Place image in the corresponding part of the combined image
        combinedImage(row*rows + (1:rows), col*cols + (1:cols), :) = resizedImage;
    end
    
    % Write the combined image to file
    imwrite(combinedImage, outputFile);
    
    fprintf('Combined image saved as %s\n', outputFile);
end



% function do_combine_images(folderPath, outputFile)
%     % Combine images from a specified folder into one single image with auto-resizing and a white background.
%     
%     % Parameters:
%     % folderPath - String, path to the folder containing images.
%     % outputFile - String, path and filename of the output image.
%     
%     % Read all image files from the folder
%     files = dir(fullfile(folderPath, '*.png'));  % Adjust the extension if needed
%     if isempty(files)
%         error('No images found in the specified folder.');
%     end
%     
%     % Determine number of images
%     numImages = length(files);
%     numPerRow = ceil(sqrt(numImages));  % Number of images per row and column
%     
%     % Read the first image to get dimensions
%     firstImage = imread(fullfile(folderPath, files(1).name));
%     [rows, cols, channels] = size(firstImage);
%     
%     % Figure to hold the montage for combining
%     hFig = figure('visible', 'off', 'Color', 'white');
%     set(hFig, 'Position', [100, 100, cols * numPerRow, rows * numPerRow]);
%     ax = axes('Parent', hFig, 'Units', 'pixels', 'Position', [0 0 cols * numPerRow rows * numPerRow]);
%     set(ax, 'Color', 'white');
%     
%     % Place each image into the combined image
%     for k = 1:numImages
%         % Determine position in grid
%         row = floor((k-1) / numPerRow);
%         col = mod(k-1, numPerRow);
%         
%         % Read and resize image
%         thisImage = imread(fullfile(folderPath, files(k).name));
%         resizedImage = imresize(thisImage, [rows, cols]);  % Resize the image
%         
%         % Place image in the corresponding part of the combined image, adjusting YData if necessary
%         image('CData', resizedImage, 'XData', [col * cols, (col + 1) * cols], 'YData', [(row + 1) * rows, row * rows], 'Parent', ax);
%     end
%     
%     % Ensure axes fit exactly to the image
%     axis(ax, 'image');
%     axis(ax, 'off');
%     
%     % Capture the figure as an image
%     frame = getframe(ax);
%     combinedImage = frame.cdata;
%     
%     % Write the combined image to file
%     imwrite(combinedImage, outputFile);
%     
%     % Close the figure
%     close(hFig);
%     
%     fprintf('Combined image saved as %s\n', outputFile);
% end


% function do_combine_images(folderPath, outputFile)
% % Combine images from a specified folder into one single image with auto-resizing.
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
% % Determine number of images
% numImages = length(files);
% numPerRow = ceil(sqrt(numImages));  % Number of images per row and column
% 
% % Read the first image to get dimensions
% firstImage = imread(fullfile(folderPath, files(1).name));
% [rows, cols, channels] = size(firstImage);
% 
% % Initialize the combined image
% combinedImage = zeros(rows * numPerRow, cols * numPerRow, channels, 'like', firstImage);
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
% % Write the combined image to file
% imwrite(combinedImage, outputFile);
% 
% fprintf('Combined image saved as %s\n', outputFile);
% end

% function do_combine_images(folderPath, outputFile)
%     % Combine images from a specified folder into one single image with auto-resizing, labeling, and a white background.
%     
%     % Parameters:
%     % folderPath - String, path to the folder containing images.
%     % outputFile - String, path and filename of the output image.
%     
%     % Read all image files from the folder
%     files = dir(fullfile(folderPath, '*.png'));  % Adjust the extension if needed
%     if isempty(files)
%         error('No images found in the specified folder.');
%     end
%     
%     % Determine number of images
%     numImages = length(files);
%     numPerRow = ceil(sqrt(numImages));  % Number of images per row and column
%     
%     % Read the first image to get dimensions
%     firstImage = imread(fullfile(folderPath, files(1).name));
%     [rows, cols, channels] = size(firstImage);
%     
%     % Figure to hold the montage for labeling
%     hFig = figure('visible', 'off', 'Color', 'white');
%     set(hFig, 'Position', [100, 100, cols * numPerRow, rows * numPerRow]);
%     ax = axes('Parent', hFig, 'Units', 'pixels', 'Position', [0 0 cols * numPerRow rows * numPerRow]);
%     set(ax, 'Color', 'white');
%     
%     % Place each image into the combined image
%     for k = 1:numImages
%         % Determine position in grid
%         row = floor((k-1) / numPerRow);
%         col = mod(k-1, numPerRow);
%         
%         % Read and resize image
%         thisImage = imread(fullfile(folderPath, files(k).name));
%         resizedImage = imresize(thisImage, [rows, cols]);  % Resize the image
%         
%         % Compute position for text label
%         xPos = (col * cols) + 0.02 * cols; % Slightly offset from the left edge of the image
%         yPos = (row * rows) + 0.98 * rows; % Slightly offset from the bottom edge of the image
%         
%         % Place image in the corresponding part of the combined image
%         image('CData', resizedImage, 'XData', [col * cols, (col + 1) * cols], 'YData', [row * rows, (row + 1) * rows], 'Parent', ax);
%         
%         % Add text label to the image
%         text(xPos, yPos, sprintf('#%d', k), 'Color', 'black', 'VerticalAlignment', 'bottom', 'FontSize', 12, 'FontWeight', 'bold');
%     end
%     
%     % Ensure axes fit exactly to the image
%     axis(ax, 'image');
%     axis(ax, 'off');
%     
%     % Capture the figure with labels as an image
%     frame = getframe(ax);
%     labeledImage = frame.cdata;
%     
%     % Write the labeled combined image to file
%     imwrite(labeledImage, outputFile);
%     
%     % Close the figure
%     close(hFig);
%     
%     fprintf('Labeled combined image saved as %s\n', outputFile);
% end

