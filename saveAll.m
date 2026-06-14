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

function  saveAll(FolderName,save_pdf,save_png,save_eps)
if nargin < 1
    FolderName = [pwd '/'];
end
if nargin < 2
    save_pdf = true;
    save_png = true;
    save_eps = true;
end
save_transparent_png = false;

FigList = findobj(allchild(0), 'flat', 'Type', 'figure');

for iFig = 1:length(FigList)
    FigHandle = FigList(iFig);
    FigName   = strcat(FolderName,strcat('fig_',get(FigHandle, 'Name')));

    savefig(FigHandle,strcat(FigName,'.fig'));
    
    if save_pdf
        exportgraphics(FigHandle,strcat(FigName,'.pdf'), 'ContentType','image',...
            'Resolution',300, 'BackgroundColor', 'none');
    end

    if save_eps
        exportgraphics(FigHandle,strcat(FigName,'.eps'),'ContentType', 'vector', ...
            'Resolution',1200, 'BackgroundColor', 'none');
    end

    if save_png
        if save_transparent_png
            set(FigHandle, 'Color', 'none'); %#ok % Set the figure background to none
            ax = get(FigHandle, 'CurrentAxes');  % Get the current axes handle
            set(ax, 'Color', 'none'); % Set the axes background to none
            export_fig(FigHandle,strcat(FigName, '.png'), '-png', '-transparent', '-r300'); 
        else
            exportgraphics(FigHandle,strcat(FigName,'.png'),'ContentType', 'image', ...
            'Resolution',300);
        end

    end

end

removeSpacesInFilenames(FolderName);