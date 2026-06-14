%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discription:                                                              %
% This is a class for evaluating von Mises stress p-norm and computing its  %
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

classdef densityVonMisesPNormElasticity < functional

    properties (GetAccess = 'public', SetAccess = 'protected')
        m_p_vm = 6; % p-norm exponent
        m_q_vm = 0.5; % stress relaxation exponent
        
        m_adjointVariable;
        m_adjointRHS;

        m_DvmDs; % derivative of von Mises stress w.r.t. stress components

        m_T1_vm;
        m_T2_vm;
        m_dSvmdx; % sensitivity of von Mises stress w.r.t design variables
        m_dpn_dvms; % derivative of p-norm with respect to von Mises stress
        m_beta_vm;
    end

    methods
        %% CONSTRUCTOR
        function obj = densityVonMisesPNormElasticity(solver, ub)
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

        function obj = setPNormExponent(obj,p)
            obj.m_p_vm = p;
        end
        % evaluate the compliance value at design and state variables queried from the elasticity solver
        % input: obj
        % output: obj, compliance
        % ubound is optional, if not provided, the functional is considered as objective and returns fx
        % if ubound is provided, the functional is considered as constraint and returns gx := fx/ubound - 1 <= 0
        function [obj,value] = evaluate(obj)

            scenarioId = 1;

            vm = obj.m_solver.m_vonMisesElems(:,:,scenarioId); % von Mises stress
            X = obj.m_solver.m_design; % design variables (pseudo-densities)

            relaxed_vm = (X.^obj.m_q_vm) .* vm; % relaxed von Mises stress

            obj.m_fx = (sum(relaxed_vm.^obj.m_p_vm,"all"))^(1/obj.m_p_vm); % p-norm

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
            % compute the adjoint right hand side
            obj = obj.computeVonMisesPNormAdjointRHS();

            % solve the adjoint problem
            obj = obj.solveAdjoint();

            % compute the gradient
            scenarioId = 1;
            X = obj.m_solver.m_design; % design variables (pseudo-densities)
            U = obj.m_solver.m_sol; % displacement field
            obj.m_T1_vm = obj.m_dpn_dvms*obj.m_beta_vm;
            obj.m_T2_vm = zeros(size(obj.m_T1_vm));
            for ely = 1:obj.m_solver.m_ny
                for elx = 1:obj.m_solver.m_nx
                    if (~obj.m_solver.m_existingElems(ely,elx)), continue; end
                    n1 = (obj.m_solver.m_ny+1)*(elx-1)+ely;
                    n2 = (obj.m_solver.m_ny+1)* elx   +ely;
                    if obj.m_solver.m_vectorize==0
                        index = [2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2];
                    else
                        elem = (obj.m_solver.m_ny)*(elx-1)+ely;
                        index = obj.m_solver.m_edofMat(elem,:);
                    end
                    KE = obj.m_solver.m_KE;
                    x = X(ely,elx);
                    Ue = U(index,scenarioId);

                    interpCoeffGrad = obj.m_solver.getInterpolationCoefficientGrad(x);

                    obj.m_T2_vm(ely,elx)=-interpCoeffGrad * obj.m_adjointVariable(index)'*KE*Ue;
                end
            end
            obj.m_dfdx = obj.m_T1_vm + obj.m_T2_vm;

            if (obj.m_isConstraint)
                grad = obj.m_dfdx/obj.m_upperBound;
            else
                grad = obj.m_dfdx;
            end
        end
    end

    methods (Access = 'private')

        function obj = solveAdjoint(obj)
            obj.m_adjointVariable = zeros(obj.m_solver.m_numDOFs,1); % initialize adjoint variable
            obj.m_adjointVariable(obj.m_solver.m_freeDOFs) = obj.m_solver.m_K(obj.m_solver.m_freeDOFs,obj.m_solver.m_freeDOFs)\obj.m_adjointRHS(obj.m_solver.m_freeDOFs); % solve linear system
            obj.m_adjointVariable(obj.m_solver.m_fixedDOFs)= 0; % set fixed DOFs to zero
        end

        function obj = computeVonMisesPNormAdjointRHS(obj)
            % Compute derivative of von Mises stress w.r.t. stress components
            scenarioId = 1;
            X = obj.m_solver.m_design; % design variables (pseudo-densities)
            vm = obj.m_solver.m_vonMisesElems(:,:,scenarioId);

            relaxed_vm = (X.^obj.m_q_vm) .* vm;
            obj.m_dpn_dvms=(sum(relaxed_vm.^obj.m_p_vm,"all"))^(1/obj.m_p_vm-1);

            nelx = obj.m_solver.m_nx;
            nely = obj.m_solver.m_ny;
            nElements = nelx * nely;
            obj.m_DvmDs = zeros(nElements,3);

            for elx = 1:nelx
                for ely = 1:nely
                    if (~obj.m_solver.m_existingElems(ely,elx)),continue;end
                    x = X(ely,elx);
                    elem = ((elx-1)*(nely) + ely);
                    Svm = relaxed_vm(ely,elx);
                    S = (x^obj.m_q_vm) * obj.m_solver.m_stressTensor(elem,:,:,scenarioId);
                    S = squeeze(S);
                    sxx = S(1,1);
                    syy = S(2,2);
                    sxy = S(1,2);

                    obj.m_DvmDs(elem,1)=1/2/Svm*(2*sxx-syy);
                    obj.m_DvmDs(elem,2)=1/2/Svm*(2*syy-sxx);
                    obj.m_DvmDs(elem,3)=3/Svm*(sxy);
                end
            end


            obj.m_adjointRHS = zeros(obj.m_solver.m_numDOFs,1); % initialize adjoint rhs
            U = obj.m_solver.m_sol; % displacement field

            % compute B
            [xi_GQ,eta_GQ,wt_GQ]= obj.m_solver.GaussQuad();
            NCell = cell(1,length(xi_GQ));
            gradNCell = cell(1,length(xi_GQ));
            for i = 1:length(xi_GQ)
                [NCell{i},gradNCell{i}] = obj.m_solver.QuadShapeFunction(xi_GQ(i),eta_GQ(i));
            end
            xNodes = obj.m_solver.m_hx*[0,1,1,0];
            yNodes = obj.m_solver.m_hy*[0,0,1,1];
            Belem = zeros(3,8);
            for g = 1:length(wt_GQ)
                gradN = gradNCell{g};
                J = obj.m_solver.Jacobian(xNodes,yNodes,xi_GQ(g),eta_GQ(g));
                dJ = det(J);
                T = J'\gradN;
                B = zeros(3,8);
                B(1,:) = [T(1,1) 0 T(1,2) 0 T(1,3) 0 T(1,4) 0];
                B(2,:) = [0 T(2,1) 0 T(2,2) 0 T(2,3) 0 T(2,4)];
                B(3,:) = [T(2,1) T(1,1) T(2,2) T(1,2) T(2,3) T(1,3) T(2,4) T(1,4)];

                Belem = Belem + wt_GQ(g)*dJ*B;
            end

            obj.m_beta_vm = zeros(size(X));
            for ely = 1:obj.m_solver.m_ny
                for elx = 1:obj.m_solver.m_nx
                    if (~obj.m_solver.m_existingElems(ely,elx)), continue; end
                    elem = ((elx-1)*(nely) + ely);
                    n1 = (obj.m_solver.m_ny+1)*(elx-1)+ely;
                    n2 = (obj.m_solver.m_ny+1)* elx   +ely;
                    index = [2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2];
                    matId = obj.m_solver.m_materialIndices(ely,elx);
                    D = obj.m_solver.getElasticityMatrix(matId);
                    x = X(ely,elx);

                    % adjoint rhs for T2
                    obj.m_adjointRHS(index) =  obj.m_adjointRHS(index) + ...
                        x^obj.m_q_vm*obj.m_dpn_dvms*Belem'*D'*obj.m_DvmDs(elem,:)'*relaxed_vm(ely,elx)^(obj.m_p_vm-1);

                    % beta for T2
                    Ue = U([2*n1-1;2*n1; 2*n2-1;2*n2; 2*n2+1;2*n2+2; 2*n1+1;2*n1+2],scenarioId);
                    obj.m_beta_vm(ely,elx)=obj.m_q_vm*x^(obj.m_q_vm-1)*relaxed_vm(ely,elx)^(obj.m_p_vm-1)*obj.m_DvmDs(elem,:)*D*Belem*Ue;
                end
            end
        end

    end
end
