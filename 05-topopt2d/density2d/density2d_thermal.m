%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Descirption:                                                              %
% This class implements the density-based topology optimization for         %
% 2D thermal problems. It inherits from the density2d class and implements  %
% the methods to solve the adjoint equations.                               %
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

classdef density2d_thermal < density2d
    methods
        function obj = density2d_thermal(solver, ...
                objective,constraints,mfgConstraints, ...
                update, ...
                maxNumIters,exportGIF, testMode)
            % check inputs
            if ~isa(solver, 'fea2d_thermal')
                error('solver must be an instance of fea2d_thermal class!');
            end

            % set default values
            if (nargin < 5),update = 'OC';end
            if (nargin < 6),maxNumIters = 250;end
            if (nargin < 7),exportGIF = false;end
            if (nargin < 8),testMode = false;end
            % construct
            obj = obj@density2d(solver,objective,constraints,mfgConstraints, ...
                update, ...
                maxNumIters,exportGIF,testMode); % call superclass

            % history
            obj.m_history.constraint.volFrac = zeros(obj.m_maxNumIters,1);
            obj.m_history.state.temperature = zeros(obj.m_maxNumIters,obj.m_solver.m_numScenarios);
        end
        %% HISTORY
        function obj = saveHistory(obj)
            obj.m_history.constraint.volFrac(obj.m_iter) = sum(obj.m_xPhys(:))/sum(obj.m_solver.m_existingElems(:));
            obj.m_history.objective(obj.m_iter)= obj.m_fx;
            obj.m_history.change(obj.m_iter) = obj.m_change;
            for scenarioId = 1:obj.m_solver.m_numScenarios
                obj.m_history.state.temperature(obj.m_iter,scenarioId) = max(obj.m_solver.m_T(:,:,scenarioId),[],'all');
            end
        end
        %% OUTPUT
        function obj = printResults(obj)
            disp(['Iteration: ' sprintf('%4i',obj.m_iter) ...
                ', Obj.: ' sprintf('%1.2e',obj.m_history.objective(obj.m_iter)) ...
                ', V/V0: ' sprintf('%1.2f',obj.m_history.constraint.volFrac(obj.m_iter))...
                ', C/C0: ' sprintf('%1.2f',obj.m_history.objective(obj.m_iter)/obj.m_history.objective(1)) ...
                ', Design Change: ' sprintf('%1.4f',obj.m_history.change(obj.m_iter))]);

            obj.m_solver.printThermalResults();
            fprintf('\n');
        end
    end
end