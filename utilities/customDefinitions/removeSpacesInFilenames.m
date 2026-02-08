%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by:                                          %
% - Amir M. Mirzendehdel, Aerospace Engineering Department, KU              %
% - Krishnan Suresh, Mechanical Engineering Department, UW-Madison          %
%                                                                           %
% Please send your comments to: amirzend@ku.edu                             %
%                                                                           %
% The code is intended for educational purposes and theoretical details     %
% are discussed in the textbook:                                            %
% Introduction to Shape and Topology Optimization using MATLAB              %
%                                                                           %
% Disclaimer:                                                               %
% The authors reserves all rights but do not guaranty that the code is      %
% free from errors. Furthermore, we shall not be liable in any event        %
% caused by the use of the program.                                         %
%                                                                           %
% License:                                                                  %
% This software is used, copied and distributed under the licensing         %
% agreement contained in the file LICENSE in the top directory of           %
% the distribution.                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function removeSpacesInFilenames(folderPath)
    % Remove spaces from filenames in a folder recursively (ignores .fig files)
    
    % Get a list of all files and directories in the folder
    filesAndDirs = dir(folderPath);
    
    % Loop through each item in the folder
    for i = 1:length(filesAndDirs)
        % Get the current file or folder name
        currentName = filesAndDirs(i).name;
        
        % Ignore '.' and '..' directories
        if strcmp(currentName, '.') || strcmp(currentName, '..')
            continue;
        end
        
        % Construct the full path of the current file or folder
        fullPath = fullfile(folderPath, currentName);
        
        % Check if the current item is a directory
        if filesAndDirs(i).isdir
            % Recursively call the function for subdirectories
            removeSpacesInFilenames(fullPath);
        else
            % Ignore .fig files
            [~, ~, ext] = fileparts(currentName);
            if ~strcmp(ext, '.fig')
                % Remove spaces from the filename
                newName = strrep(currentName, ' ', '');
                
                % If the name changed, rename the file
                if ~strcmp(newName, currentName)
                    movefile(fullPath, fullfile(folderPath, newName));
                end
            end
        end
    end
end
