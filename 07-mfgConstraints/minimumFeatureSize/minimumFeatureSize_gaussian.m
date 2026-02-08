%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discription:                                                              %
% This is a class for imposing minimum feature size constraints             %
% based on convolution kernel used in shape and topology optimization       %
% (primarily used in level-set methods.)                                    %
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

classdef  minimumFeatureSize_gaussian < mfgConstraints
    properties (GetAccess = 'public', SetAccess = 'protected')
        m_sigma; % standard deviation
        m_filter; % gaussian filter (created through MATLAB's 'fspecial')
    end
    methods
        %% CONSTRUCTOR
        function obj = minimumFeatureSize_gaussian(solver,sigma)
            % check if solver is valid
            if (~isa(solver, 'simulation2d')), error('solver must be an instance of simulation2d class!');end % check if solver is valid

            % constructor based on superclass
            obj = obj@mfgConstraints(solver);

            if (nargin > 2)
                obj.m_sigma = sigma;
            else
                obj.m_sigma = 0.6;
            end
            % feature size filter
            obj.m_filter = fspecial('gaussian', [3 3],obj.m_sigma); % smoothen topological sensitivity field
        end


        function [filteredDesign] = filterDesign(obj, design)
            % filter the design variables
            % input: obj, design variables
            % output: obj, filtered design variables
            filteredDesign = design .* obj.m_solver.m_existingElems;
        end

        function [filteredSensitivity] = filterSensitivity(obj, ~, sensField,nIter)
            % filter the sensitivity fields
            % input: obj, sensitivity fields
            % output: obj, filtered sensitivity fields
            if nargin < 4, nIter = 2; end
            for i=1:nIter
                filteredSensitivity = filter2(obj.m_filter,sensField); % smoothen the field
                filteredSensitivity = filteredSensitivity .* obj.m_solver.m_existingElems; % zero sensitivity on inactive regions
            end
        end
    end
end
