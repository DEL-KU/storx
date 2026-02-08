%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discription:                                                              %
% This is a class for imposing edge retain constraints used in              %
% shape and topology optimization (e.g., SIMP).                             %
%                                                                           %
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

classdef  symmetry_density < mfgConstraints
    properties (GetAccess = 'public', SetAccess = 'protected')
        m_dir; % 0: x-dir, 1: y-dir
    end
    methods
        %% CONSTRUCTOR
        function obj = symmetry_density(solver,dir)
            % check if solver is valid
            if (~isa(solver, 'simulation2d')), error('solver must be an instance of simulation2d class!');end % check if solver is valid

            % constructor based on superclass
            obj = obj@mfgConstraints(solver);


            obj.m_dir = dir;
            if (obj.m_dir~=0 && obj.m_dir~=1)
                error('symmetry direction must be 0 (x-dir) or 1 (y-dir)!');
            end
        end

        function [filteredDesign] = filterDesign(obj, design)
            % filter the design variables
            % input: obj, design variables
            % output: obj, filtered design variables

            if obj.m_dir == 0 % x-dir
                filteredDesign = max(design , fliplr(design));
            else % y-dir
                filteredDesign = max(design , flipud(design));
            end
        end

        function [filteredSensitivity] = filterSensitivity(obj, ~, sensField)
            % filter the sensitivity fields
            % input: obj, sensitivity fields
            % output: obj, filtered sensitivity fields

            if obj.m_dir == 0 % x-dir
                filteredSensitivity = min(sensField , fliplr(sensField));
            else % y-dir
                filteredSensitivity = min(sensField , flipud(sensField));
            end
        end
    end
end
