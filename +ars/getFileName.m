%%% extract file name from path
function fileName = getFileName(file)
    [~,name,ext] = fileparts(file);
    fileName = [name ext];
end