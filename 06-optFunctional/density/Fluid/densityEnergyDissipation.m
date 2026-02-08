%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discription:                                                              %
% This is a class for evaluating fluid energy dissipation and computing its               %
% gradient for density-based topology optimization.                         %
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

classdef densityEnergyDissipation < functional
    properties (GetAccess = 'public', SetAccess = 'protected')
        m_adjointRHS; % right hand side of the adjoint problem
        m_adjointVariable; % adjoints for the fluid problem
    end

    methods
        %% CONSTRUCTOR
        function obj = densityEnergyDissipation(solver, ub)
            % check if solver is valid
            if (~isa(solver, 'fea2d_fluid')), error('solver must be an instance of fea2d_fluid class!');end

            % if upper bound is provided set the value, otherwise set it to NaN
            upper_bound = NaN;
            if (nargin > 1)
                upper_bound = ub;
            end

            % constructor based on superclass functional
            obj = obj@functional(solver, upper_bound);
        end

        % evaluate the compliance value at design and state variables queried from the elasticity solver
        % input: obj
        % output: obj, compliance
        % ubound is optional, if not provided, the functional is considered as objective and returns fx
        % if ubound is provided, the functional is considered as constraint and returns gx := fx/ubound - 1 <= 0
        function [obj,value] = evaluate(obj)
            obj.m_fx = zeros(obj.m_solver.m_numScenarios,1);
            dxv = obj.m_solver.m_hx*ones(1,obj.m_solver.m_numElems);
            dyv = obj.m_solver.m_hy*ones(1,obj.m_solver.m_numElems);
            muv = obj.m_solver.m_materials(1).mu*ones(1,obj.m_solver.m_numElems); % Assuming single fluid
            for scenarioId = 1:obj.m_solver.m_numScenarios
                U = obj.m_solver.m_sol(:,scenarioId);
                obj.m_fx(scenarioId) = sum(PHI(dxv,dyv,muv, ...
                    obj.m_solver.m_alpha(:)', ...
                    U(obj.m_solver.m_edofMat')));
            end

            if (obj.m_numEvaluations == 0)
                obj.m_fx0 = obj.m_fx;
                obj.m_scale = obj.m_fx0;
            end

            if (obj.m_isConstraint)
                value = obj.m_fx/obj.m_upperBound - 1; % as constraint
            else
                value = obj.m_fx; % as objective
            end
            obj.m_numEvaluations = obj.m_numEvaluations + 1; % update number of evaluations
        end

        % gradient of the compliance at design and state variables queried from the solver.
        % additional adjoint problems are solved in this function to obtain the lagrnage multipliers
        % input: obj, ubound (optional)
        % output: obj, functional gradient w.r.t. design variables
        % ubound is optional, if not provided, the functional is considered as objective and returns df/dx
        % if ubound is provided, the functional is considered as constraint and returns dgx/dx := dfx/dx/ubound
        function [obj,grad] = gradient(obj)
            % solve adjoint problem
            obj = obj.solveAdjoint();

            % compute gradient
            obj.m_dfdx = zeros(obj.m_solver.m_ny,obj.m_solver.m_nx,obj.m_solver.m_numScenarios);
            dxv = obj.m_solver.m_hx*ones(1,obj.m_solver.m_numElems);
            dyv = obj.m_solver.m_hy*ones(1,obj.m_solver.m_numElems);
            muv = obj.m_solver.m_materials(1).mu*ones(1,obj.m_solver.m_numElems); % Assuming single fluid
            rhov = obj.m_solver.m_materials(1).rho*ones(1,obj.m_solver.m_numElems); % Assuming single fluid

            obj.m_solver = obj.m_solver.computeInterpolationCoefficientGrad(obj.m_solver.m_design);

            for scenarioId = 1:obj.m_solver.m_numScenarios
                U = obj.m_solver.m_sol(:,scenarioId);
                sR = dRESdg(dxv,dyv,muv,rhov, ...
                    obj.m_solver.m_alpha(:)', ...
                    obj.m_solver.m_alphaGrad(:)', ...
                    U(obj.m_solver.m_edofMat'));

                dRdg = sparse(obj.m_solver.m_iR(:),obj.m_solver.m_jE(:),sR(:));

                dphidg = dPHIdg(dxv,dyv,muv, ...
                    obj.m_solver.m_alpha(:)', ...
                    obj.m_solver.m_alphaGrad(:)', ...
                    U(obj.m_solver.m_edofMat'));

                obj.m_dfdx(:,:,scenarioId) = reshape(dphidg - obj.m_adjointVariable'*dRdg, ...
                    obj.m_solver.m_ny,obj.m_solver.m_nx);
            end

            if (obj.m_isConstraint)
                grad = obj.m_dfdx/obj.m_upperBound;
            else
                grad = obj.m_dfdx;
            end
        end
    end

    methods(Access = 'private')
        function obj = solveAdjoint(obj)

            obj.m_solver = obj.m_solver.assembleK();

            obj.m_adjointRHS = zeros(obj.m_solver.m_numDOFs,1);
            obj.m_adjointVariable = zeros(obj.m_solver.m_numDOFs,1);

            dxv = obj.m_solver.m_hx*ones(1,obj.m_solver.m_numElems);
            dyv = obj.m_solver.m_hy*ones(1,obj.m_solver.m_numElems);
            muv = obj.m_solver.m_materials(1).mu*ones(1,obj.m_solver.m_numElems); % Assuming single fluid

            U = obj.m_solver.m_sol;
            sR = [dPHIds(dxv,dyv,muv,obj.m_solver.m_alpha(:)',U(obj.m_solver.m_edofMat')); zeros(4,obj.m_solver.m_numElems)];
            obj.m_adjointRHS = sparse(obj.m_solver.m_iR,obj.m_solver.m_jR,sR(:));
            obj.m_adjointRHS(obj.m_solver.m_fixedDOFs) = 0;
            obj.m_adjointVariable  = obj.m_solver.m_K'\obj.m_adjointRHS;
        end
    end
end
