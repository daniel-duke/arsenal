%%% recursively remove folder
function myRmdir(dirPath)
    if ~isfolder(dirPath)
        error('The specified path is not a directory or does not exist.');
    end
    
    files = dir(dirPath);
    for i = 1:length(files)
        if strcmp(files(i).name, '.') || strcmp(files(i).name, '..')
            continue;
        end
        filePath = fullfile(dirPath, files(i).name);
        if isfolder(filePath)
            myRmdir(filePath);
        else
            delete(filePath);
        end
    end
    rmdir(dirPath);
end
