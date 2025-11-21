%%% search given folder for references to arsenal functions
function ars_refs_allFiles = findArsReferences(searchFold)

    %%% set default search folder as current folder
    if nargin == 0
        searchFold = "./";
    end

    %%% location of arsenal folder
    arsFold = "/Users/dduke/Programs/arsenal/+ars/";

    %%% get a list of all .m files in the arsenal folder
    arsFiles = dir(fullfile(arsFold, '*.m'));
    
    %%% extract the names of all available arsenal functions
    ars_functions = cell(length(arsFiles), 1);
    for k = 1:length(arsFiles)
        [~, searchFileName, ~] = fileparts(arsFiles(k).name);
        ars_functions{k} = searchFileName;
    end
    ars_functions = unique(ars_functions);

    %%% get a list of all .m files in the search folder (recursively)
    searchFiles = dir(fullfile(searchFold, '**', '*.m'));

    %%% initialize the results structure
    ars_refs_allFiles = struct();

    %%% loop through each file in the search folder
    for k = 1:length(searchFiles)
        [~, searchFileName, ~] = fileparts(searchFiles(k).name);

        %%% get full file path
        searchFile = fullfile(searchFiles(k).folder, searchFiles(k).name);
        
        %%% read the content of the file
        content = fileread(searchFile);
        
        %%% get a list of all function references in the file
        all_refs_oneFile = extractRefs(content);

        %%% filter only the arsenal functions
        ars_refs_oneFile = intersect(all_refs_oneFile, ars_functions);

        %%% remove arsenal functions defining themselves
        removed = 0;
        for i = 1:length(ars_refs_oneFile)
            if strcmp(ars_refs_oneFile{i-removed},searchFileName)
                ars_refs_oneFile(i-removed) = [];
                removed = removed + 1;
            end
        end

        %%% store the result if there are any arsenal references
        if ~isempty(ars_refs_oneFile)
            ars_refs_allFiles.(searchFileName) = ars_refs_oneFile;
        end
    end

    %%% display result
    fprintf("\n")
    files = fieldnames(ars_refs_allFiles);
    for i = 1:length(files)
        disp(files{i})
        references = ars_refs_allFiles.(files{i});
        for j = 1:length(references)
            fprintf("- %s\n",references{j})
        end
    end
    fprintf("\n")
end

function references = extractRefs(content)
    %%% match function calls using a regular expression
    pattern = '(?<![\w])([a-zA-Z][a-zA-Z0-9_]*)\s*(?=\()';
    matches = regexp(content, pattern, 'match');
    
    %%% clean up matches to get function names only
    references = unique(strtrim(regexprep(matches, '\s*\(', '')));
    
    %%% exclude MATLAB built-in keywords
    matlabKeywords = {
        'if', 'else', 'elseif', 'end', 'for', 'while', 'switch', 'case', ...
        'otherwise', 'try', 'catch', 'function', 'return', 'break', ...
        'continue', 'global', 'persistent', 'parfor', 'spmd', 'classdef', ...
        'methods', 'properties', 'events'};
    references = setdiff(references, matlabKeywords);
end
