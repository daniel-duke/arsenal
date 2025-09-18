%%% create folder, deleting old folder if it exists
function createEmptyFold(dirPath)
    if isfolder(dirPath)
        ars.myRmdir(dirPath);
    end
    mkdir(dirPath)
end