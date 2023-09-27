function checkOrCreateDir(directory)
    if ~exist(directory, 'dir')
        mkdir(directory);
        disp('Folder created successfully.');
    else
        disp('Folder already exists.');
    end
end
