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
function  export_gifs(FolderName)
    % EXPORT_GIFS Exports all open figures as GIFs to a specified folder
    %
    % Input Arguments:
    %     FolderName - (optional) directory where GIFs will be saved
    %                  defaults to the current working directory if not provided

    if nargin < 1, FolderName = [pwd '/']; end  % Set default folder if none provided
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');  % Find all open figure handles
    for iFig = 1:length(FigList)  % Loop through each figure
        FigHandle = FigList(iFig);  % Get the current figure handle
        FigName   = strcat(FolderName,strcat('fig_',get(FigHandle, 'Name')));  % Construct the filename
        ax = get(FigHandle, 'CurrentAxes');  % Get the current axes handle
        title(ax, '');  % Remove the title by setting it to an empty string
        exportgraphics(FigHandle,strcat(FigName,'.gif'),'Append',true);  % Export figure as GIF
    end

