function releaseLock(lockFilePath)
    if exist(lockFilePath, 'file')
        delete(lockFilePath);
    end
end
