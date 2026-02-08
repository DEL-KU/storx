%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                              %
% This class implements a 2D pareto-tracing topology optimization for       %
% elasticity problems using topological sensitivity fields and              %
% fixed-point iteration. It inherits from the pareto2d class.               %
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

classdef pareto2d_elasticity < pareto2d

    methods
        function obj = pareto2d_elasticity(solver,objective,constraints, mfgConstraints, ...
                volDecrement,paretoAggressiveness,exportGIF,testMode)

            % set default values
            if (nargin < 5),volDecrement = 0.025;end
            if (nargin < 6),paretoAggressiveness = 0.65;end
            if (nargin < 7),exportGIF = false;end
            if (nargin < 8),testMode = false;end
            % construct
            obj = obj@pareto2d(solver,objective,constraints, mfgConstraints, ...
                volDecrement,paretoAggressiveness,exportGIF,testMode) % call superclass
        end

        %% OBJECTIVE & CONSTRAINTS
        function obj = evaluate(obj)
            % evaluate objective
            [obj.m_objective,obj.m_fx] = obj.m_objective.evaluate();
            obj.m_fx = mean(obj.m_fx); % take mean if needed for number of scenarios
            if (obj.m_iter==1), obj.m_fx0 = obj.m_fx;    end

            % evaluate constraints
            obj.m_gx = zeros(obj.m_numConstraints,1);
            for g = 1 : obj.m_numConstraints
                [obj.m_constraints{g}, gx] = obj.m_constraints{g}.evaluate();
                obj.m_gx(g) = gx;
            end
        end
        %% GRADIENT OVERRIDE
        function obj = gradient(obj)
            % gradient objective
            obj.m_dfdx = zeros(obj.m_solver.m_ny,obj.m_solver.m_nx,obj.m_solver.m_numScenarios);
            [obj.m_objective,obj.m_dfdx] = obj.m_objective.gradient();
            obj.m_dfdx = mean(obj.m_dfdx,3);
            % gradient constraints
            obj.m_dgdx = [];
            for g = 1 : obj.m_numConstraints
                [obj.m_constraints{g}, dgdx] = obj.m_constraints{g}.gradient();
                dgdx = dgdx(:); % convert to vector
                obj.m_dgdx = [obj.m_dgdx;dgdx']; % append
            end
        end

        %% HISTORY
        function obj = saveHistory(obj)
            obj.m_history.constraint.volFrac(obj.m_iter) = sum(obj.m_x(:))/sum(obj.m_solver.m_existingElems(:));
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
                ', C/C0: ' sprintf('%1.2f',obj.m_history.objective(obj.m_iter)/obj.m_history.objective(1)) ...
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

