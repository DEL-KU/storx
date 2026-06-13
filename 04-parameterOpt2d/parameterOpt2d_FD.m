%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                              %
% Finite-difference parameter optimization (concrete subclass).             %
% For fea2d_elasticity solvers: uses fmincon with internal FD gradients.   %
% For triFEA2d_elasticity solvers: uses a semi-analytic gradient computed  %
% via mesh morphing (snapNodesToBRep) and the adjoint sensitivity           %
%   C' = f^T u',  K u' = -K' u.                                            %
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
                obj.m_totalEvals = 30;
            else
                obj.m_totalEvals = 30 * (1 + obj.m_numParams);
            end
            obj.openProgressBar('FD: Finite Difference');
            obj.m_feasibleExploredSolutions = struct('X', [], 'Fval', []);

            LB = obj.m_param0.lb ./ obj.m_param0.value;
            UB = obj.m_param0.ub ./ obj.m_param0.value;
            x0 = ones(size(LB));

            if isa(obj.m_solverInitial, 'triFEA2d_elasticity')
                % Semi-analytic gradient via mesh morphing
                opt = optimoptions('fmincon', 'Display', 'iter', ...
                    'TolX',                     obj.m_terminationTolerance, ...
                    'TolFun',                   obj.m_terminationTolerance, ...
                    'ConstraintTolerance',      obj.m_terminationTolerance, ...
                    'SpecifyObjectiveGradient', true);

                [xMin, ~, flag, output] = fmincon( ...
                    @obj.evaluateObjectiveWithGradient, x0, [], [], [], [], ...
                    LB, UB, @obj.evaluateConstraint, opt);
            else
                % Direct FD: fmincon approximates the gradient internally
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
            % Semi-analytic (indirect) FD objective + gradient for triFEA2d_elasticity.
            % Uses mesh projection via snapNodesToBRep to morph the baseline mesh to
            % the perturbed geometry, then solves K u' = -K' u for the displacement
            % sensitivity, giving C' = f^T u'.

            params = x .* obj.m_param0.value;
            da     = obj.m_finiteDifferenceStepSize;

            fem0 = obj.solveQuiet(params);
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
            f0      = fem0.m_F;
            u0      = fem0.m_Sol;
            freeDOF = fem0.m_FreeDOF;
            K0_free = K0(freeDOF, freeDOF);
            nDOF    = length(u0);

            p0_backup    = fem0.m_mesh.p;
            brep0_backup = fem0.m_brep;

            grad = zeros(size(x));
            for i = 1:obj.m_numParams
                da_phys        = da * obj.m_param0.value(i);
                params_pert    = params;
                params_pert(i) = params(i) + da_phys;
                brep_pert      = obj.m_brepHandle(params_pert);

                fem0.m_brep = brep_pert;
                fem0 = fem0.snapNodesToBRep();
                fem0 = fem0.assembleK();
                K_pert = fem0.m_K;

                Kprime = (K_pert - K0) / da_phys;

                rhs          = -Kprime(freeDOF, :) * u0;
                u_prime_free = K0_free \ rhs;

                u_prime          = zeros(nDOF, 1);
                u_prime(freeDOF) = u_prime_free;

                C_prime = f0' * u_prime;
                grad(i) = C_prime * obj.m_param0.value(i) / obj.m_cx0;

                fem0.m_mesh.p = p0_backup;
                fem0.m_brep   = brep0_backup;
            end

            fem0.m_K = K0;
            fem0.m_F = f0;
        end

    end
end
