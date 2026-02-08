%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decription:                                                               %
% This class Implements a 2D level-set topology optimization for elasticity %
% problems. It inherits from the standardHJ2d class and incorporates        %
% topological senstivity in the Hamilton-Jacobi function.                   %
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

classdef modifiedHJ2d_elasticity < standardHJ2d
    properties
        m_dTdx; % topological sensitivity of compliance function
        m_topWeight; % weigth of topological sensitivity
    end

    methods
        function obj = modifiedHJ2d_elasticity(solver,objective,constraints, mfgConstraints, ...
                topWeight, ...
                maxNumIters,exportGIF,testMode)

            % set default values
            if (nargin < 5),topWeight = 10;end
            if (nargin < 6),maxNumIters = 300;end
            if (nargin < 7),exportGIF = false;end
            if (nargin < 8),testMode = false;end

            numReinit = 5;
            stepLength = 30;
            nHolesX = 0;   nHolesY = 0;   r0 = 0;
            % construct
            obj = obj@standardHJ2d(solver,objective,constraints, mfgConstraints, ...
                nHolesX,nHolesY,r0, ...
                numReinit,stepLength,maxNumIters,exportGIF,testMode);

            % Topological Sensitivity Weight
            obj.m_topWeight = topWeight;

            % history
            obj.m_history.constraint.volFrac = zeros(obj.m_maxNumIters,1);
            obj.m_history.state.deformation = zeros(obj.m_maxNumIters,obj.m_solver.m_numScenarios);
            obj.m_history.state.vonMises = zeros(obj.m_maxNumIters,obj.m_solver.m_numScenarios);
        end

        %% UPDATE DESIGN
        function obj = update(obj)
            ve = min(obj.m_solver.m_ve,1e-3);
            % Update Lagrange multipliers for augmented Lagrangian
            obj = obj.updateLagrangeMultipliers();

            % Combine sensitivities based on augmented Lagrangian
            obj.m_dSdx = obj.m_dfdx.dSdx - ve *obj.m_lagLower + ...
                ve *1/obj.m_lagUpper*obj.m_gx;
            obj.m_dTdx = obj.m_dfdx.dTdx + ve *pi*(obj.m_lagLower ...
                - 1/obj.m_lagUpper*obj.m_gx);

            % Smooth/filter the sensitivities
            obj = obj.filterSensitivity();

            % Ensure only the values for existing elements are used
            obj.m_dSdx = obj.m_solver.m_existingElems.*obj.m_dSdx;
            obj.m_dTdx = obj.m_solver.m_existingElems.*obj.m_dTdx;

            % Extend sensitivites using a zero border
            vFull = zeros(size(obj.m_dSdx)+2); vFull(2:end-1,2:end-1) = -obj.m_dSdx;
            gFull = zeros(size(obj.m_dTdx)+2); gFull(2:end-1,2:end-1) = obj.m_dTdx;

            % Choose time step for evolution based on CFL value
            dt = 0.1/(obj.m_stepLength*max(abs(obj.m_dSdx(:))));

            % Evolve for total time stepLength * CFL value:
            for i = 1:obj.m_stepLength
                % Calculate derivatives on the grid
                dpx = circshift(obj.m_lsf,[0,-1])-obj.m_lsf;
                dmx = obj.m_lsf - circshift(obj.m_lsf,[0,1]);
                dpy = circshift(obj.m_lsf,[-1,0]) - obj.m_lsf;
                dmy = obj.m_lsf - circshift(obj.m_lsf,[1,0]);
                % Update level set function using an upwind scheme
                obj.m_lsf = obj.m_lsf - dt * min(vFull,0).* ...
                    sqrt( min(dmx,0).^2+max(dpx,0).^2+min(dmy,0).^2+max(dpy,0).^2 ) ...
                    - dt * max(vFull,0) .*...
                    sqrt( max(dmx,0).^2+min(dpx,0).^2+max(dmy,0).^2+min(dpy,0).^2 )...
                    - obj.m_topWeight*dt*gFull;
            end
            % Save the current design
            xOld = obj.m_x;

            % New structure obtained from lsf
            obj.m_lsf(2:end-1,2:end-1) = obj.m_solver.m_existingElems.*obj.m_lsf(2:end-1,2:end-1);
            lsf = obj.m_lsf(2:end-1,2:end-1);
            obj.m_x = obj.m_solver.m_existingElems.*(lsf<0);

            % Filter the density
            obj = obj.filterDensity();

            % evaluate change
            delta_x = abs(xOld(obj.m_solver.m_existingElems==1) - obj.m_x(obj.m_solver.m_existingElems==1));
            obj.m_change = sum(delta_x(:))/obj.m_solver.m_numExistingElems;
        end
        %% FILTERS
        function obj = filterDensity(obj)
            % filter density
            for mfgConsId = 1:obj.m_numMfgConstraints
                obj.m_x = obj.m_mfgConstraints{mfgConsId}.filterDesign(obj.m_x);
            end
        end

        function obj = filterSensitivity(obj)
            % filter sensitivity of objective and constraints
            % the inputs are assumed to be flattened colummn vectors
            % each number of the vector corresponds to an element of the mesh
            % each nx by ny matrix corresponds to a scenario that is reshpaed to a matrix
            % of size nx by ny and then filtered and flattened again
            % The retain filter is applied to the each sensitivity field here
            % as well as the ensuring that we only have values for existing elements

            for mfgConsId = 1:obj.m_numMfgConstraints
                % Augment Shape Sensitivity
                obj.m_dSdx = obj.m_mfgConstraints{mfgConsId}.filterSensitivity(obj.m_x,obj.m_dSdx);
                obj.m_dTdx = obj.m_mfgConstraints{mfgConsId}.filterSensitivity(obj.m_x,obj.m_dTdx);
                % Apply existing elements
                obj.m_dSdx = obj.m_dSdx .* obj.m_solver.m_existingElems;
                obj.m_dTdx = obj.m_dTdx .* obj.m_solver.m_existingElems;
            end
        end
        
        %% OBJECTIVE & CONSTRAINTS
        function obj = evaluate(obj)
            % evaluate objective
            if (isa(obj.m_objective,'modifiedHJComplianceElasticity'))
                [obj.m_objective,obj.m_fx] = obj.m_objective.evaluate();
                obj.m_fx = mean(obj.m_fx);
            else
                disp('Only modifiedHJComplianceElasticity objective is available for this optimization!');
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
            if (isa(obj.m_objective,'modifiedHJComplianceElasticity'))
                % m_dfdx has shape (m_dfdx.dSdx) and topology (m_dfdx.dTdx) sensitivities of objective
                [obj.m_objective,obj.m_dfdx] = obj.m_objective.gradient();
            else
                disp('Only modifiedHJComplianceElasticity objective is available for this optimization!');
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

