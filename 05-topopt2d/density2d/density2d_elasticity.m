%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                              %
% This class implements the density-based topology optimization (e.g., SIMP)%
% for 2D elasticity problems. It inherits from the density2d class and      %
% implements the methods for evaluating and computing the gradient of       %
% the compliance,von Mises stress p-norm, and displacement objectives.      %
% It also handles the  adjoint variable and sensitivity                     %
% calculations for these objectives.                                        %
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

classdef density2d_elasticity < density2d
    methods
        function obj = density2d_elasticity(solver,...
                objective,constraints,mfgConstraints, ...
                update, ...
                maxNumIters,exportGIF,testMode)

            % check inputs
            if ~isa(solver, 'fea2d_elasticity')
                error('solver must be an instance of fea2d_elasticity class');
            end
            % set default values
            if (nargin < 4)
                mfgConstraints.rmin = 1.5;
            end
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
            obj.m_history.state.deformation = zeros(obj.m_maxNumIters,obj.m_solver.m_numScenarios);
            obj.m_history.state.vonMises = zeros(obj.m_maxNumIters,obj.m_solver.m_numScenarios);
        end
        %% HISTORY
        function obj = saveHistory(obj)
            obj.m_history.constraint.volFrac(obj.m_iter) = sum(obj.m_xPhys(:))/sum(obj.m_solver.m_existingElems(:));
            obj.m_history.objective(obj.m_iter)= obj.m_fx;
            obj.m_history.change(obj.m_iter) = obj.m_change;
            for scenarioId = 1:obj.m_solver.m_numScenarios
                obj.m_history.state.deformation(obj.m_iter,scenarioId) = max(obj.m_solver.m_def(:,:,scenarioId),[],'all');
                obj.m_history.state.vonMises(obj.m_iter,scenarioId) = max(obj.m_solver.m_vonMisesElems(:,:,scenarioId),[],'all');
            end
        end
        %% OUTPUT
        function obj = printResults(obj)
            disp(['Iteration: ' sprintf('%4i',obj.m_iter) ...
                ', Obj.: ' sprintf('%1.2e',obj.m_history.objective(obj.m_iter)) ...
                ', V/V0: ' sprintf('%1.2f',obj.m_history.constraint.volFrac(obj.m_iter))...
                ', f/f0: ' sprintf('%1.2f',obj.m_history.objective(obj.m_iter)/obj.m_history.objective(1)) ...
                ', Design Change: ' sprintf('%1.4f',obj.m_history.change(obj.m_iter))]);

            for scenarioId = 1:obj.m_solver.m_numScenarios
                disp(['Scenario: ' sprintf('%4i',scenarioId) ...
                    ', Deformation: ' sprintf('%1.2e',obj.m_history.state.deformation(obj.m_iter,scenarioId)) ...
                    ', vonMises: ' sprintf('%1.2e',obj.m_history.state.vonMises(obj.m_iter,scenarioId)) ...
                    ])
            end
            fprintf('\n');
        end
    end
end