%%% determine if file exists
function found = checkFileExist(file,name,required)
    arguments
        file char
        name char = 'the'
        required = true
    end
    found = true;
    f = fopen(file,'r');
    if f == -1
        if required == true
            error('Could not find %s file:\n%s',name,file);
        else
            sprintf('Could not find %s file:\n%s',name,file);
            found = false;
        end
    end
end