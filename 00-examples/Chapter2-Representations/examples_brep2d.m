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

clc; clear;  close all;format compact; format long

brepClass = @brep2d;
example = 1;
switch (example)
    case 1
        %% Beam with Chamfers
        geom = brep2d('BeamWithChamfers.brep');
    case 2
        %% Beam with Fillet
        BeamWithFillets.vertices = [0 0;1.7 0;2 0.3;2 0.7;1.7 1;0 1;1.7 0.3;1.7 0.7]';
        BeamWithFillets.segments = [1 1 2 0;2 2 3 -7;1 3 4 0;2 4 5 -8;1 5 6 0;1 6 1 0]';
        geom = brep2d(BeamWithFillets);
    case 3
        %% Beam with Fillet
        BeamWithHole.vertices = [0 0;1 0;1 0.3;1 0.7;2 0;2 1;0 1;1 0.5]';
        BeamWithHole.segments = [1 1 2 0; -1 2 3 0; 2 3 4 8; 2 4 3 8; -1 3 2 0; 1 2 5 0; 1 5 6 0; 1 6 7 0; 1 7 1 0]';
        geom = brep2d(BeamWithHole);
end

geom.plotGeometryWithLabels();