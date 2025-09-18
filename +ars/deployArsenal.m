%%% update functions in deployed arsenal directory
function deployArsenal(arsDepFold)

    %%% set default deployed arsenal folder
    if nargin == 0
        arsDepFold = "./+ars/";
    end

    %%% location of arsenal folder
    arsFold = "/Users/dduke/Files/analysis/arsenal/+ars/";

    %%% get a list of all arsenal files
    arsFiles = dir(arsFold);

    %%% iterate through each arsenal file
    for i = 1:length(arsFiles)
        arsFile = arsFiles(i).name;
        if arsFiles(i).isdir || arsFile(1) == '.'
            continue;
        end

        %%% copy, only if file exists in deployed arsenal
        if exist(fullfile(arsDepFold, arsFile), 'file') && arsFile ~= "deployArsenal.m"
            copyfile(fullfile(arsFold, arsFile), fullfile(arsDepFold, arsFile));
        end
    end
end
