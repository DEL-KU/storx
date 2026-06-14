%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                              %
% Finite-difference parameter optimization (concrete subclass).             %
% For fea2d_elasticity solvers: uses fmincon with internal FD gradients.   %
% For triFEA2d_elasticity solvers: uses a semi-analytic gradient computed  %
% via mesh projection. With constant loads, the compliance derivative is    %
%   C' = -u^T K' u.                                                        %
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

classdef parameterOpt2d_FD < parameterOpt2d

    methods

        function obj = parameterOpt2d_FD(brepHandle, solverHandle, param, ...
                objective, constraints, ...
                terminationTolerance, finiteDifferenceStepSize, ...
                exportGIF, testMode)

            if nargin < 8 || isempty(exportGIF), exportGIF = false; end
            if nargin < 9 || isempty(testMode),  testMode  = false; end

            obj = obj@parameterOpt2d(brepHandle, solverHandle, param, ...
                objective, constraints, ...
                terminationTolerance, finiteDifferenceStepSize, ...
                exportGIF, testMode);
        end

        function obj = optimize(obj)
            if isa(obj.m_solverInitial, 'triFEA2d_elasticity')
                obj.m_totalEvals = 100;
            else
                obj.m_totalEvals = 100 * (1 + obj.m_numParams);
            end
            obj.openProgressBar('FD: Finite Difference');
            obj.m_feasibleExploredSolutions = struct('X', [], 'Fval', []);

            LB = obj.m_param0.lb ./ obj.m_param0.value;
            UB = obj.m_param0.ub ./ obj.m_param0.value;
            x0 = ones(size(LB));

            if isa(obj.m_solverInitial, 'triFEA2d_elasticity')
                scenarioId = 1;
                obj.m_solverInitial = obj.solveProjectedTri(obj.m_param0.value);
                obj.m_cx0        = obj.m_solverInitial.computeCompliance();
                obj.m_maxDef0    = obj.m_solverInitial.m_maxDef(scenarioId);
                obj.m_maxStress0 = obj.m_solverInitial.m_maxStress(scenarioId);

                % Semi-analytic gradient with SQP for small parameter spaces.
                opt = optimoptions('fmincon', ...
                'Display', 'iter', ...
                'Algorithm', 'sqp', ...
                'StepTolerance',       obj.m_terminationTolerance, ...
                'FunctionTolerance',   obj.m_terminationTolerance, ...
                'OptimalityTolerance', obj.m_terminationTolerance, ...
                'ConstraintTolerance', obj.m_terminationTolerance, ...
                'SpecifyObjectiveGradient',  true, ...
                'SpecifyConstraintGradient', false, ...
                'FiniteDifferenceType',      'central', ...
                'FiniteDifferenceStepSize',  obj.m_finiteDifferenceStepSize);

                [xMin, ~, flag, output] = fmincon( ...
                    @obj.evaluateObjectiveWithGradient, x0, [], [], [], [], ...
                    LB, UB, @obj.evaluateConstraint, opt);
            else
                % Direct FD: fmincon approximates gradients internally.
                opt = optimoptions('fmincon', 'Display', 'iter', ...
                    'TolX',                     obj.m_terminationTolerance, ...
                    'TolFun',                   obj.m_terminationTolerance, ...
                    'ConstraintTolerance',      obj.m_terminationTolerance, ...
                    'FiniteDifferenceStepSize', obj.m_finiteDifferenceStepSize);

                [xMin, ~, flag, output] = fmincon(@obj.evaluateObjective, x0, [], [], ...
                    [], [], LB, UB, @obj.evaluateConstraint, opt);
            end

            % fmincon reports (1 + nParams) FEA calls per iteration for FD gradient
            funcCount = output.funcCount * (1 + obj.m_numParams);
            obj = obj.finalizeResults(xMin, flag, funcCount);

            if ~obj.m_testMode
                obj = obj.plotConvergence();
            end
        end

        function [fx, grad] = evaluateObjectiveWithGradient(obj, x)
            % Semi-analytic objective gradient for triFEA2d_elasticity.
            % Mesh projection keeps the connectivity fixed while K' is formed.

            params = x .* obj.m_param0.value;

            fem0 = obj.solveProjectedTri(params);
            cx   = fem0.computeCompliance();
            fx   = cx / obj.m_cx0;

            obj.m_paramsHistory    = [obj.m_paramsHistory;    params];
            obj.m_objectiveHistory = [obj.m_objectiveHistory; cx];
            inBounds = all(params <= obj.m_param0.ub & params >= obj.m_param0.lb);
            obj.m_feasibleHistory  = [obj.m_feasibleHistory;  double(inBounds)];
            if ~obj.m_testMode
                n = numel(obj.m_objectiveHistory);
                fprintf('  eval %-4d | C = %.6f\n', n, cx);
            end
            obj.updateProgressBar(cx);
            if obj.m_exportGIF, export_gifs(); end

            % Extract baseline FEA quantities
            K0      = fem0.m_K;
            u0      = fem0.m_Sol;
            freeDOF = fem0.m_FreeDOF;
            u0_free = u0(freeDOF);

            p0_backup    = fem0.m_mesh.p;
            brep0_backup = fem0.m_brep;

            grad = zeros(size(x));
            for i = 1:obj.m_numParams
                [h, direction] = obj.boundSafePhysicalStep(params, i);
                if h == 0
                    grad(i) = 0;
                    continue;
                end

                params_step = params;
                params_step(i) = params(i) + direction * h;
                K_step = obj.projectedStiffness(fem0, params_step, p0_backup, brep0_backup);

                if direction > 0
                    Kprime = (K_step - K0) / h;
                else
                    Kprime = (K0 - K_step) / h;
                end

                C_prime = -u0_free' * Kprime(freeDOF, freeDOF) * u0_free;
                grad(i) = C_prime * obj.m_param0.value(i) / obj.m_cx0;
            end

            fem0.m_K = K0;
            fem0.m_mesh.p = p0_backup;
            fem0.m_brep   = brep0_backup;
        end

        function [h, direction] = boundSafePhysicalStep(obj, params, paramId)
            h = obj.m_finiteDifferenceStepSize * abs(obj.m_param0.value(paramId));
            if h <= 0
                h = sqrt(eps) * max(1, abs(params(paramId)));
            end

            forwardRoom  = obj.m_param0.ub(paramId) - params(paramId);
            backwardRoom = params(paramId) - obj.m_param0.lb(paramId);

            if forwardRoom >= h
                direction = 1;
            elseif backwardRoom >= h
                direction = -1;
            elseif forwardRoom > 0
                h = forwardRoom;
                direction = 1;
            elseif backwardRoom > 0
                h = backwardRoom;
                direction = -1;
            else
                h = 0;
                direction = 0;
            end
        end

        function K = projectedStiffness(obj, fem, params, p0, brep0)
            fem.m_mesh.p = p0;
            fem.m_brep   = obj.m_brepHandle(params);
            fem          = obj.projectMeshToBRep(fem, p0);
            fem          = fem.assembleK();
            K            = fem.m_K;
            fem.m_mesh.p = p0;
            fem.m_brep   = brep0;
        end

        function fem = solveProjectedTri(obj, params)
            fem = obj.createSolver(params);
            fem = obj.projectMeshToBRep(fem, fem.m_mesh.p);
            fem = fem.assembleK();
            fem = fem.assembleBC();
            fem = fem.solve();
            fem = fem.postProcess();
        end

        function fem = projectMeshToBRep(~, fem, p0)
            bndryNodes = unique(fem.m_mesh.e(1:fem.m_nodesPerEdge, :));
            bndryPts   = p0(:, bndryNodes);
            [~, projectedPts] = fem.distOfPointsToBrep(bndryPts);

            bndryDisp = projectedPts - bndryPts;
            p = p0;
            if numel(bndryNodes) >= 3
                Fx = scatteredInterpolant(bndryPts(1,:)', bndryPts(2,:)', ...
                    bndryDisp(1,:)', 'natural', 'nearest');
                Fy = scatteredInterpolant(bndryPts(1,:)', bndryPts(2,:)', ...
                    bndryDisp(2,:)', 'natural', 'nearest');
                p = p0 + [Fx(p0(1,:)', p0(2,:)')'; Fy(p0(1,:)', p0(2,:)')'];
            end
            p(:, bndryNodes) = projectedPts;
            fem.m_mesh.p = p;
        end

    end

end
