function do_combine_images(folderPath, outputFile, layout)
% do_combine_images  Combine images into a single file with flexible layout.
%
% Usage:
%   do_combine_images(folderPath, outputFile)          % default square grid
%   do_combine_images(folderPath, outputFile, 'row')   % one horizontal row
%   do_combine_images(folderPath, outputFile, 'column')% one vertical column
%
% Inputs
%   folderPath - string, path to folder with .png images
%   outputFile - string, name of the merged image to create
%   layout     - (optional) 'grid' | 'row' | 'column'   (default: 'grid')
%
% The script lets you reorder images interactively, then stitches them
% together on a white background using the specified layout.

    if nargin < 3 || isempty(layout),  layout = 'grid';  end
    layout = validatestring(lower(layout), {'grid','row','column'}, ...
                            mfilename, 'layout', 3);

    % ---------------------------------------------------------------------
    % 1. Gather image files
    files = dir(fullfile(folderPath, '*.png'));
    if isempty(files)
        error('No PNG images found in "%s".', folderPath);
    end

    [~, idx] = sort({files.name});
    files = files(idx);
    files = files(~strcmp({files.name}, outputFile));

    % Optional re-order step (unchanged):
    fprintf('Images found:\n');
    for i = 1:numel(files)
        fprintf('%2d: %s\n', i, files(i).name);
    end
    inp = input( ...
        '\nEnter comma-separated order (or press Enter to keep): ','s');
    if ~isempty(inp)
        newOrd = str2num(inp); %#ok<ST2NM>
%         if numel(unique(newOrd)) ~= numel(files)
%             error('Invalid order.');
%         end
        files = files(newOrd);
    end

    % ---------------------------------------------------------------------
    % 2. Determine canvas size based on layout
    nImg = numel(files);
    ref   = imread(fullfile(folderPath, files(1).name));
    [h, w, c] = size(ref);

    switch layout
        case 'row'
            rows = 1;     cols = nImg;
        case 'column'
            rows = nImg;  cols = 1;
        otherwise  % 'grid'
            cols = ceil(sqrt(nImg));
            rows = ceil(nImg / cols);
    end

    blank = 255 * ones(h*rows, w*cols, c, class(ref));  % white background

    % ---------------------------------------------------------------------
    % 3. Paste each image
    for k = 1:nImg
        img = imresize(imread(fullfile(folderPath, files(k).name)), [h w]);

        r0 = floor((k-1)/cols)*h + 1;
        c0 = mod(k-1, cols)*w + 1;
        blank(r0:r0+h-1, c0:c0+w-1, :) = img;
    end

    % ---------------------------------------------------------------------
    % 4. Save
    outPath = fullfile(folderPath, outputFile);
    imwrite(blank, outPath);
    fprintf('Merged image saved to %s\n', outPath);
end

% function do_combine_images(folderPath, outputFile)
% % do_combine_images Combines images from a specified folder into one single image 
% % with auto-resizing, a white background, and optional manual re-ordering.
% %
% % Usage:
% %   do_combine_images(folderPath, outputFile)
% %
% % Inputs:
% %   folderPath - (String) Path to the folder containing images (e.g., PNGs).
% %   outputFile - (String) Name of the final merged file (e.g., 'combined.png').
% %
% % The function lists all images found, gives you an opportunity to
% % reorder them (by typing a custom sequence of indices), and merges
% % them into a grid saved in outputFile.
% 
%     % 1. Read image files
%     files = dir(fullfile(folderPath, '*.png'));  % Adjust extension if needed
%     if isempty(files)
%         error('No images found in the specified folder.');
%     end
% 
%     % 2. Sort files by name (alphabetically), then exclude outputFile
%     [~, idx] = sort({files.name});
%     files = files(idx);
%     files = files(~strcmp({files.name}, outputFile));
% 
%     % 3. Sort remaining files by modification date (ascending)
%     fileDates = [files.datenum];
%     [~, idx] = sort(fileDates);
%     files = files(idx);
%     
%     % 4. Display the files to be merged
%     fprintf('The following images will be merged (current order):\n');
%     for i = 1:length(files)
%         fprintf('%2d: %s\n', i, files(i).name);
%     end
%     
%     % 5. Ask user for a custom re-ordering
% %     promptMsg = sprintf( ...
% %         '\nEnter a comma-separated list of indices (e.g. 3,1,2,...) to re-order,\n' ...
% %         'or press Enter to keep the current order: ');
% %     userInput = input(promptMsg, 's');
%     
%     promptMsg = sprintf('\nEnter a comma-separated list of indices (e.g. 3,1,2,...) to re-order,\nor press Enter to keep the current order: ');
%     
%     userInput = input(promptMsg, 's');
%     
%     
%     if ~isempty(userInput)
%         % Parse user input into numeric array
%         newOrder = str2num(userInput); %#ok<ST2NM> 
%         if isempty(newOrder) || length(newOrder) ~= length(files) ...
%                 || any(newOrder < 1) || any(newOrder > length(files)) ...
%                 || length(unique(newOrder)) ~= length(files)
%             error(['Invalid indices. Ensure you provide a comma-separated list ' ...
%                    'with each index from 1 to %d exactly once.'], length(files));
%         end
%         files = files(newOrder);
%         fprintf('Images re-ordered by user input.\n');
%     else
%         fprintf('Proceeding with default order.\n');
%     end
% 
%     % 6. Print final order
%     fprintf('\nFinal merge order:\n');
%     for i = 1:length(files)
%         fprintf('%2d: %s\n', i, files(i).name);
%     end
%     fprintf('\n');
% 
%     % 7. Merge images into a grid
%     numImages = length(files);
%     numPerRow = ceil(sqrt(numImages));  % # of images per row/column
% 
%     % Read first image to get reference dimensions
%     firstImage = imread(fullfile(folderPath, files(1).name));
%     [rows, cols, channels] = size(firstImage);
% 
%     % Initialize white background
%     if isa(firstImage, 'uint8')
%         combinedImage = 255 * ones(rows*numPerRow, cols*numPerRow, channels, 'uint8');
%     else
%         combinedImage = 255 * ones(rows*numPerRow, cols*numPerRow, channels, 'like', firstImage);
%     end
% 
%     % 8. Place each image in the combined grid
%     for k = 1:numImages
%         % Grid position (row/col in the final big image)
%         gridRow = floor((k-1) / numPerRow);
%         gridCol = mod(k-1, numPerRow);
% 
%         % Read and resize the image to match reference dimensions
%         thisImage = imread(fullfile(folderPath, files(k).name));
%         resizedImage = imresize(thisImage, [rows, cols]);
% 
%         % Insert it into combinedImage
%         rowStart = gridRow*rows + 1;
%         colStart = gridCol*cols + 1;
%         combinedImage(rowStart:(rowStart+rows-1), ...
%                       colStart:(colStart+cols-1), :) = resizedImage;
%     end
% 
%     % 9. Write the combined image
%     export_path = fullfile(folderPath, outputFile);
%     imwrite(combinedImage, export_path);
% 
%     % 10. Confirm
%     fprintf('Combined image saved as: %s\n', export_path);
%     cd(folderPath);
% end
% 
