function isLocked = acquireLock(lockFilePath)
    if exist(lockFilePath, 'file')
        % Lock file exists, do not proceed with analysis
        isLocked = false;
    else
        % No lock file, attempt to create one to proceed
        fclose(fopen(lockFilePath, 'w'));
        isLocked = true;
    end
end
