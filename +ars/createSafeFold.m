%%% create folder, only if it doesn't exist
function createSafeFold(dirPath)
    if ~isfolder(dirPath)
        mkdir(dirPath)
    end
end