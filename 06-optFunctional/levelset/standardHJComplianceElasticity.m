%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discription:                                                              %
% This is a class for evaluating compliance and computing its               % 
% gradient for level-set shape optimization using standard Hamilton-Jacobi  %
% equation.                                                                 %
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

classdef standardHJComplianceElasticity < functional
    properties (GetAccess = 'public', SetAccess = 'protected')
        m_ce; % compliance at every element
    end

    methods
        %% CONSTRUCTOR
        function obj = standardHJComplianceElasticity(solver, ub)
            % check if solver is valid
            if (~isa(solver, 'fea2d_elasticity')), error('solver must be an instance of fea2d_elasticity class!');end

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
            KE = obj.m_solver.m_KE;
            
            if (obj.m_solver.m_vectorize==0)
                U = obj.m_solver.m_sol;
                for scenarioId = 1:obj.m_solver.m_numScenarios
                    for ely = 1:obj.m_solver.m_ny
                        for elx = 1:obj.m_solver.m_nx
                            if (~obj.m_solver.m_existingElems(ely,elx)), continue; end
                            n1 = (obj.m_solver.m_ny+1)*(elx-1)+ely;
                            n2 = (obj.m_solver.m_ny+1)* elx   +ely;
                            Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],scenarioId);
                            obj.m_fx(scenarioId) = obj.m_fx(scenarioId) + Ue'*KE*Ue;
                        end
                    end
                end
            else
                obj.m_ce = zeros(obj.m_solver.m_ny,obj.m_solver.m_nx,obj.m_solver.m_numScenarios);
                X = obj.m_solver.m_design;
                for scenarioId = 1:obj.m_solver.m_numScenarios
                    U = obj.m_solver.m_sol(:,scenarioId);
                    % obj.m_ce(:,:,scenarioId) = reshape(sum((U(obj.m_solver.m_edofMat)*KE).*U(obj.m_solver.m_edofMat),2),obj.m_solver.m_ny,obj.m_solver.m_nx);
                    % obj.m_fx(scenarioId) = sum(sum(X.*obj.m_ce(:,:,scenarioId)));
                    F = obj.m_solver.m_f(:,scenarioId);
                    obj.m_fx(scenarioId) = F'*U;
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
            % Set material parameters, find Lame values
            dSdx = zeros(obj.m_solver.m_ny,obj.m_solver.m_nx,obj.m_solver.m_numScenarios);
            U = obj.m_solver.m_sol;
            X = obj.m_solver.m_design;
            for scenarioId = 1:obj.m_solver.m_numScenarios
                for ely = 1:obj.m_solver.m_ny
                    for elx = 1:obj.m_solver.m_nx
                        if (~obj.m_solver.m_existingElems(ely,elx)), continue; end
                        elem = (obj.m_solver.m_ny)*(elx-1)+ely;
                        n1 = (obj.m_solver.m_ny+1)*(elx-1)+ely;
                        n2 = (obj.m_solver.m_ny+1)* elx   +ely;
                        if (obj.m_solver.m_vectorize==0)
                            Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],scenarioId);
                        else
                            Ue = U(obj.m_solver.m_edofMat(elem,:),scenarioId);
                        end
                        dSdx(ely,elx,scenarioId) = -max(X(ely,elx),0.001)*Ue'*obj.m_solver.m_KE*Ue;
                    end
                end
            end
            dSdx = mean(dSdx,3); % take average of all scenarios
            dSdx = dSdx / sum(obj.m_fx0);

            obj.m_dfdx = dSdx; % shape sensitivity

            if (obj.m_isConstraint)
                grad = obj.m_dfdx/obj.m_upperBound;
            else
                grad = obj.m_dfdx;
            end
        end
    end
end
