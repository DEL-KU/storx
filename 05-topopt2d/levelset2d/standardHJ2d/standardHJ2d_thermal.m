%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                              %
% This class implements a 2D level-set based topology optimization using    %
% the standard Hamilton-Jacobi (HJ) method for thermal problems.            %    
%                                                                           %
% The level-set approach and code is based on Challis, Vivien J.            %
% "A discrete level-set topology optimization code written in Matlab."      %
% Structural and multidisciplinary optimization 41 (2010): 453-464.         %
%                                                                           %
% For higher order solution, the LSF re-initialization to SDF is achieved   %
% using the "ToolboxLS" developed by:                                       %
% Ian M. Mitchell (mitchell@cs.ubc.ca), 2004                                %
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

classdef standardHJ2d_thermal < standardHJ2d

    methods
        function obj = standardHJ2d_thermal(solver,objective,constraints, mfgConstraints, ...
                nHolesX,nHolesY,r0, ...
                maxNumIters,exportGIF,testMode)

            % set default values
            % Hole initializtion
            if (nargin < 5),nHolesX = 0;end
            if (nargin < 6),nHolesY = 0;end
            if (nargin < 7),r0 = 0;end
            if (nargin < 8),maxNumIters = 300;end
            if (nargin < 9),exportGIF = false;end
            if (nargin < 10),testMode = false;end


            numReinit = 5;
            stepLength = 30;
            % construct
            obj = obj@standardHJ2d(solver,objective,constraints, mfgConstraints, ...
                nHolesX,nHolesY,r0, ...
                numReinit,stepLength,maxNumIters,exportGIF,testMode); % call superclass

             % history
            obj.m_history.constraint.volFrac = zeros(obj.m_maxNumIters,1);
            obj.m_history.state.temperature = zeros(obj.m_maxNumIters,obj.m_solver.m_numScenarios);
        end
        %% OBJECTIVE & CONSTRAINTS
        function obj = evaluate(obj)
            % evaluate objective
            if (isa(obj.m_objective,'standardHJComplianceThermal'))
                [obj.m_objective,obj.m_fx] = obj.m_objective.evaluate();
                obj.m_fx = mean(obj.m_fx);
            else
                disp('Only standardHJComplianceThermal objective is available for this optimization!');
            end
            if (obj.m_iter == 1), obj.m_fx0 = obj.m_fx; end

            % evaluate constraints
            obj.m_gx = 0; % only volume constraint is implemented
            for g = 1 : obj.m_numConstraints
                if (isa(obj.m_constraints{g},'volume'))
                    [obj.m_constraints{g}, gx] = obj.m_constraints{g}.evaluate();
                    obj.m_gx = gx*obj.m_constraints{g}.m_upperBound;
                    break;
                else
                    disp('Only volume constraint is available for this optimization!');
                end
            end
        end

        function obj = gradient(obj)
            % evaluate objective
            if (isa(obj.m_objective,'standardHJComplianceThermal'))
                [obj.m_objective,obj.m_dfdx] = obj.m_objective.gradient();
            else
                disp('Only standardHJComplianceThermal objective is available for this optimization!');
            end

            % evaluate constraints
            obj.m_dgdx = zeros(size(obj.m_x));
            for g = 1 : obj.m_numConstraints
                if (isa(obj.m_constraints{g},'volume'))
                    [obj.m_constraints{g}, dgdx] = obj.m_constraints{g}.gradient();
                    dgdx = dgdx*obj.m_constraints{g}.m_upperBound;
                    obj.m_dgdx = dgdx;
                    break;
                else
                    disp('Only volume constraint is available for this optimization!');
                end
            end
        end

        %% HISTORY
        function obj = saveHistory(obj)
            obj.m_history.constraint.volFrac(obj.m_iter) = sum(obj.m_x(:))/sum(obj.m_solver.m_existingElems(:));
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
            
            for scenarioId = 1:obj.m_solver.m_numScenarios
                disp(['Scenario: ' sprintf('%4i',scenarioId) ...
                    ', Temperature: ' sprintf('%1.2e',obj.m_history.state.temperature(obj.m_iter,scenarioId)) ...
                    ])
            end
            fprintf('\n');
        end
    end
end

