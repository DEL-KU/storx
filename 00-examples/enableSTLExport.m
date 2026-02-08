% Recursive function to search for all .m files in subfolders and process them
function enableSTLExport()
clc
% Start the recursive search in the current folder
currentScript = mfilename('fullpath');
[path,~,~] = fileparts(currentScript);
disp(path)
processFolder(path,currentScript);
disp('STL export is truned on for all examples.');
end

function processFolder(folderPath,currentScriptPath)
% Get all .m files in the current folder
files = dir(fullfile(folderPath, '*.m'));

% Process each .m file in this folder
for k = 1:length(files)
    filePath = fullfile(files(k).folder, files(k).name);

    % Skip the current script
    if strcmp(files(k).name, 'disableSTLExport.m') || strcmp(files(k).name, 'enableSTLExport.m')
        fprintf('Skipping script: %s\n', filePath);
        continue;
    end

    % Read file
    fid = fopen(filePath, 'r');
    if fid == -1
        fprintf('Could not open file: %s\n', filePath);
        return
    end
    fileContents = fread(fid, '*char')';
    fclose(fid);

    % Split into lines
    lines = regexp(fileContents, '\r\n|\n|\r', 'split');

    pattern = '\<exportSTL\>\s*=\s*(?:false|0)\s*;?';  % word-boundary safe in MATLAB
    replacement = 'exportSTL = true;';

    changed = false;
    for i = 1:numel(lines)
        L = lines{i};

        % Separate code from any comment starting with %
        pct = strfind(L, '%');
        if ~isempty(pct)
            code = L(1:pct(1)-1);
            comment = L(pct(1):end);
        else
            code = L; comment = '';
        end

        newCode = regexprep(code, pattern, replacement);
        if ~strcmp(newCode, code)
            changed = true;
        end

        lines{i} = [newCode comment];
    end

    if changed
        updatedContents = strjoin(lines, newline);
        fid = fopen(filePath, 'w');
        if fid == -1
            fprintf('Could not write to file: %s\n', filePath);
        else
            fwrite(fid, updatedContents, '*char');
            fclose(fid);
            fprintf('Processed file: %s - exportSTL enabled\n', filePath);
        end
    else
        fprintf('Processed file: %s - exportSTL already enabled or not found\n', filePath);
    end


end

% Recursively process all subfolders
subfolders = dir(folderPath);
for i = 1:length(subfolders)
    if subfolders(i).isdir && ~ismember(subfolders(i).name, {'.', '..'})
        processFolder(fullfile(folderPath, subfolders(i).name),currentScriptPath);
    end
end
end
