%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discription:                                                              %
% This is an abstract class for manufacturing constraints used in shape and %
% topology optimization.                                                    %
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

classdef (Abstract)  mfgConstraints < handle
    properties (GetAccess = 'public', SetAccess = 'protected')
        m_solver; % solver
    end
    methods (Abstract)
        % abstract methods
        % filter the design variables
        [filteredDesign] = filterDesign(obj, x)

        % filter the sensitivity fields
        [filteredSensitivity] = filterSensitivity(obj, sensField);
    end
    methods
        %% CONSTRUCTOR
        function obj = mfgConstraints(solver)
            % constructor
            if (~isa(solver, 'simulation2d')), error('solver must be an instance of simulation2d class!');end % check if solver is valid
            obj.m_solver = solver; % assign solver
        end
    end
end
