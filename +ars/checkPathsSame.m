%%% determine if two paths point to the same location
function isSame = checkPathsSame(path1,path2)
    path1 = expandTilde(path1);
    path2 = expandTilde(path2);

    path1 = getAbsolutePath(path1);
    path2 = getAbsolutePath(path2);

    path1 = java.io.File(path1);
    path2 = java.io.File(path2);
    isSame = strcmp(char(path1.getCanonicalPath()), char(path2.getCanonicalPath()));
end

function expanded = expandTilde(path)
    if startsWith(path, '~')
        if ispc
            home = getenv('USERPROFILE');
        else
            home = getenv('HOME');
        end
        expanded = strrep(path, '~', home);
    else
        expanded = path;
    end
end

function absPath = getAbsolutePath(path)
    if startsWith(path, './') || startsWith(path, '../') || ~startsWith(path, '/')
        absPath = fullfile(pwd, path);
    else
        absPath = path;
    end
end