%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discription:                                                              %
% This is a class for evaluating compliance and computing its               %
% gradient using topological sensitivity fieldsfor evolutionary and         %
% Pareto-tracing topology optimization.                                     %
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

classdef topologicalSensitivityComplianceElasticity < functional
    properties (GetAccess = 'public', SetAccess = 'protected')
        m_ce; % compliance at every element
        m_df_min;
        m_df_max;
    end

    methods
        %% CONSTRUCTOR
        function obj = topologicalSensitivityComplianceElasticity(solver, ub)
            % check if solver is valid
            if (~isa(solver, 'fea2d_elasticity')), error('solver must be an instance of fea2d_elasticity class!');end

            % if upper bound is provided set the value, otherwise set it to NaN
            upper_bound = NaN;
            if (nargin > 1)
                upper_bound = ub;
            end

            % constructor based on superclass functional
            obj = obj@functional(solver, upper_bound);
            obj.m_dfdx = zeros(obj.m_solver.m_ny,obj.m_solver.m_nx,obj.m_solver.m_numScenarios); % initialize to 0 only at the beginning
        end

        % evaluate the compliance value at design and state variables queried from the elasticity solver
        % input: obj
        % output: obj, compliance
        % ubound is optional, if not provided, the functional is considered as objective and returns fx
        % if ubound is provided, the functional is considered as constraint and returns gx := fx/ubound - 1 <= 0
        function [obj,value] = evaluate(obj)
            obj.m_fx = zeros(obj.m_solver.m_numScenarios,1);
            KE = obj.m_solver.m_KE;
            X = obj.m_solver.m_design;
            if (obj.m_solver.m_vectorize==0)
                U = obj.m_solver.m_sol;
                for scenarioId = 1:obj.m_solver.m_numScenarios
                    for ely = 1:obj.m_solver.m_ny
                        for elx = 1:obj.m_solver.m_nx
                            if (~obj.m_solver.m_existingElems(ely,elx)), continue; end
                            if (~X(ely,elx)), continue; end
                            n1 = (obj.m_solver.m_ny+1)*(elx-1)+ely;
                            n2 = (obj.m_solver.m_ny+1)* elx   +ely;
                            Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],scenarioId);
                            x = X(ely,elx);
                            interpCoeff = obj.m_solver.getInterpolationCoefficient(x);
                            KElem = interpCoeff * KE;

                            obj.m_fx(scenarioId) = obj.m_fx(scenarioId) + Ue'*KElem*Ue;
                        end
                    end
                end
            else
                xPhys = obj.m_solver.getInterpolationCoefficient(X);
                obj.m_ce = zeros(obj.m_solver.m_ny,obj.m_solver.m_nx,obj.m_solver.m_numScenarios);
                for scenarioId = 1:obj.m_solver.m_numScenarios
                    U = obj.m_solver.m_sol(:,scenarioId);
                    obj.m_ce(:,:,scenarioId) = reshape(sum((U(obj.m_solver.m_edofMat)*KE).*U(obj.m_solver.m_edofMat),2),obj.m_solver.m_ny,obj.m_solver.m_nx);
                    obj.m_fx(scenarioId) = sum(sum(xPhys.*obj.m_ce(:,:,scenarioId)));
                end
            end
            if (obj.m_numEvaluations == 0), obj.m_fx0 = obj.m_fx; end

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
            % Compute the topological sensitivity at the center of each element
            X = obj.m_solver.m_design;
            nu = obj.m_solver.m_materials(1).nu;
            for scenarioId = 1:obj.m_solver.m_numScenarios
                for  elx = 1:obj.m_solver.m_nx
                    for  ely = 1:obj.m_solver.m_ny
                        if (~obj.m_solver.m_existingElems(ely,elx)),continue;end
                        if (~X(ely,elx)), continue; end
                        elem = ((elx-1)*(obj.m_solver.m_ny) + ely);
                        stressTensor = reshape(obj.m_solver.m_stressTensor(elem,:,:,scenarioId),2,2);
                        strainTensor = reshape(obj.m_solver.m_strainTensor(elem,:,:,scenarioId),2,2);
                        obj.m_dfdx(ely,elx,scenarioId) = 4/(1+nu)*sum(sum(stressTensor.*strainTensor))-  ...
                            (1-3*nu)/(1-nu^2)*trace(stressTensor)*trace(strainTensor);
                    end
                end
            end

            if (obj.m_numEvaluations == 1)
                obj.m_df_min = min(abs(obj.m_dfdx(:))); 
                obj.m_df_max = max(abs(obj.m_dfdx(:))); 
            end

            if (obj.m_isConstraint)
                grad = obj.m_dfdx/obj.m_upperBound/obj.m_fx0;
            else
                grad = (obj.m_dfdx-obj.m_df_min)/(obj.m_df_max-obj.m_df_min);
            end
        end
    end
end
