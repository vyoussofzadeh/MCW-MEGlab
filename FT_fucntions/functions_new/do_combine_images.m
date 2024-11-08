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

% Filter out the outputFile if it exists in the same folder
files = files(~strcmp({files.name}, outputFile));

% Extract date info and sort files by date
fileDates = [files.datenum];
[~, idx] = sort(fileDates);
files = files(idx);

% Assuming you've already sorted `files` as needed
fprintf('The following images will be merged:\n');
for i = 1:length(files)
    fprintf('%d: %s\n', i, files(i).name);
end

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
imwrite(combinedImage, fullfile(folderPath,outputFile));

fprintf('Combined image saved as %s\n', outputFile);
cd(folderPath)
end
