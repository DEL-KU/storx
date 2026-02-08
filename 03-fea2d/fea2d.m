%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discription:                                                              %
% This is an abstract class for 2D finite element analysis (FEA) simulations%
% in the context of shape and topology optimization. It inherits from the   %
% simulation2d class and provides methods for pre-processing, solving, and  %
% post-processing the simulation. The class is designed to handle design    %
% variables, materials, and load scenarios in a 2D grid mesh environment.   %
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

classdef (Abstract) fea2d < simulation2d
    properties(GetAccess = 'public', SetAccess = 'protected')
        m_BCtype; % for each segment, 0 means Neumann, 1 means Dirichlet
        m_BCvalue; % for each segment, the value of Dirichlet or Neumann BC
        m_numDOFperNode; % number of degrees of freedom per node
        m_numDOFs;  % total number of degrees of freedom
        m_numDOFperElem;   % number of degrees of freedom per element
        m_freeDOFs;  % free DOFs (not fixed)
        m_fixedDOFs; % fixed DOFs
        m_forcedNodes;  % nodes with non-zero force
        m_fixedNodes; % nodes with Dirichlet BC
        m_K; % global stiffness matrix
        m_KE; % elemental stiffness matrix
        m_f; % global force vector (nDOFx1xnScenarios)
        m_fBody; % global body force
        m_fE; % elemental body force
        m_fixed; % dirichlet values
        m_sol; % solution vector
        m_penalty; % current penalty factor
        m_penaltyStruct; % range of penalty values and increment e.g., struct('min',2,'max',3,'inc',0.05);
        m_nonExistingDOF; % excluding non-existing DOF from solve
        %% Vectorized parameters
        m_edofVec; % elemental degrees of freedom vector
        m_edofMat;   % elemental degrees of freedom matrix
        m_iK;    % row indices for global stiffness matrix
        m_jK;    % column indices for global stiffness matrix
    end
    methods (Abstract)
        obj = computeMaterialPropertiesMatrices(obj)
        obj = assembleBC(obj)
        obj = assembleK(obj)
        obj = postProcess(obj)
        obj = assembleInternalLoad(obj)
    end
    methods
        function obj = fea2d(brep,numElements,numDOFperNode,materials,...
                interpolation,numScenarios,penaltyStruct,uniformGrid)
            % Constructor for the 2D finite element analysis class
            % Inputs:
            %   brep: boundary representation of the geometry
            %   numElements: number of elements in the mesh
            %   numDOFperNode: number of degrees of freedom per node (e.g., 2 for 2D elasticity)
            %   materials: material properties for the simulation
            %   interpolation: interpolation method for the finite element analysis needed for design optimization (e.g., SIMP)
            %   numScenarios: number of load scenarios to consider
            %   penaltyStruct: structure containing penalty parameters (min, max, inc)
            %   uniformGrid: boolean indicating whether to use a uniform grid (default is false)

            % construct
            obj = obj@simulation2d(brep,numElements,materials,...
                interpolation,numScenarios,uniformGrid); % call superclass

            obj.m_numDOFperNode = numDOFperNode; % e.g., 2 for 2D elasticity
            obj.m_numDOFs = obj.m_numDOFperNode*(obj.m_nx+1)*(obj.m_ny+1); % total number of degrees of freedom
            obj.m_sol = zeros(obj.m_numDOFs,numScenarios); % initialize solution
            obj.m_f = zeros(obj.m_numDOFs,numScenarios);     % initialize force vector
            obj.m_BCtype = zeros(obj.m_numBndrySegs,obj.m_numDOFperNode,obj.m_numScenarios);% e.g., [u v] for elasticity
            obj.m_BCvalue = zeros(obj.m_numBndrySegs,obj.m_numDOFperNode,obj.m_numScenarios);% e.g., [u v] for elasticity
            obj.m_penalty = penaltyStruct.min;  % initial penalty factor
            obj.m_penaltyStruct = penaltyStruct; % range of penalty values

            obj.m_nonExistingDOF = []; % excluding non-existing DOF from solve
            for nodex = 1:obj.m_nx+1
                for nodey = 1:obj.m_ny+1
                    node = (nodex-1)*(obj.m_ny+1) + nodey;
                    if obj.m_existingNodes(nodey,nodex), continue;end
                    obj.m_nonExistingDOF = [obj.m_nonExistingDOF; 2*node-1;2*node];
                end
            end
            obj.m_nonExistingDOF = unique(obj.m_nonExistingDOF);
        end
        %% PRE_PROCESS
        function obj = preProcess(obj)
            % Pre-process the finite element analysis simulation
            obj = obj.computeMaterialPropertiesMatrices();
            obj = obj.assembleBC();
        end
        %% SOLVE LINEAR SYSTEM
        function obj = solve(obj)
            % Solve the linear system of equations for the finite element analysis
            % This method assembles the global stiffness matrix, applies boundary conditions,
            % and solves for the nodal displacements.
            % It assembles the global stiffness matrix (K)
            % and assumes that the force vector (f) and boundary conditions (BC) have been initialized
            % and that the boundary conditions have been applied.
            % The solution is stored in the m_sol property.
            % The method also handles Dirichlet boundary conditions by modifying the force vector
            % and the stiffness matrix accordingly.
            % The solution is computed using a direct solver (backslash operator).
            obj = obj.assembleInternalLoad();
            obj = obj.assembleK();
            obj.m_sol = zeros(obj.m_numDOFs,obj.m_numScenarios);

            fTilde = obj.m_f; % create a local copy
            fTilde = fTilde + obj.m_fBody; % add body force
            % now subtract all the dirichlet values from rhs
            % assuming dirichlet bc is common for all scenarios
            fTilde = fTilde(obj.m_freeDOFs,:);
            for dof = 1:obj.m_numDOFs
                if (abs(obj.m_fixed(dof)) > 0)
                    fTilde = fTilde - obj.m_K(obj.m_freeDOFs,dof)*obj.m_fixed(dof);  % subtract dirichlet value
                end
            end
            % direct solve
            obj.m_sol(obj.m_freeDOFs,:) = obj.m_K(obj.m_freeDOFs,obj.m_freeDOFs)\fTilde;
            for scenarioId = 1:obj.m_numScenarios
                obj.m_sol(obj.m_fixedDOFs,scenarioId)= obj.m_fixed(obj.m_fixedDOFs); % add dirichlet value
            end
        end
        %% BOUNDARY CONDITIONS
        function obj = fixEdge(obj,boundaryEdges,value)
            % Fix the specified boundary edges with a Dirichlet condition
            % value is the value of the Dirichlet condition, default is 0
            % Inputs:
            %   boundaryEdges: vector of boundary edges to fix
            %   value: value of the Dirichlet condition, default is 0
            % Output:
            %   obj: updated object

            if nargin < 2, error('Boundary edges must be specified'); end
            if nargin < 3, value = 0; end
            obj.m_BCtype(boundaryEdges,1:obj.m_numDOFperNode,1) = 1;
            obj.m_BCvalue(boundaryEdges,1,1) = value;
        end

        function obj = applyDirichletOnDOF(obj,dof)
            if (nargin < 2), error('DOF must be specified'); end

            if (dof < 1 || dof > obj.m_numDOFs)
                error('DOF must be between 1 and %d',obj.m_numDOFs);
            end

            node = ceil(dof/obj.m_numDOFperNode); % node number
            obj.m_fixedDOF = unique([obj.m_fixedDOF dof]); % dof number
            obj.m_fixedNodes = unique([obj.m_fixedNodes node]); % node number
        end
        %% SET NODAL FORCES
        function obj = setForce(obj,f)
            % Set the nodal force vector for the finite element analysis
            % Inputs:
            %   f: force vector of size (numDOFs x numScenarios)
            % Output:
            %   obj: updated object with the force vector set
            % This method sets the nodal force vector for the finite element analysis simulation.
            % It checks the size and validity of the force vector, ensuring it matches the number of
            % degrees of freedom and scenarios defined in the simulation.
            % It also checks for NaN, Inf, and negative values in the force vector.
            % If the force vector is valid, it updates the m_f property of the object.

            if (nargin < 2), error('Force vector must be specified'); end
            if (size(f,1) ~= obj.m_numDOFs || size(f,2) ~= obj.m_numScenarios)
                error('Force vector must be of size %d x %d',obj.m_numDOFs,obj.m_numScenarios);
            end
            if (any(isnan(f(:))))
                error('Force vector contains NaN values');
            end
            if (any(isinf(f(:))))
                error('Force vector contains Inf values');
            end

            obj.m_f = f;
        end
        %% GET NODAL FORCES
        function f = getForce(obj)
            f = obj.m_f ;
        end
        %% SET NODAL SOLUTION
        function obj = setSolution(obj,sol)
            obj.m_sol = sol;
        end
        %% GET NODAL SOLUTION
        function sol = getSolution(obj)
            sol = obj.m_sol;
        end
        %% SET GLOBAL STIFFNESS MATRIX
        function obj = setK(obj,K)
            obj.m_K = K;
        end
        %% GET GLOBAL STIFFNESS MATRIX
        function K = getK(obj)
            K = obj.m_K;
        end
        %% SET PENALTY FACTOR FOR PROPERTY INTERPOLATION
        function obj = setPenaltyFactor(obj,penalty)
            obj.m_penalty = penalty;
        end
        %% GET PENALTY FACTOR FOR PROPERTY INTERPOLATION
        function penalty = getPenaltyFactor(obj)
            penalty = obj.m_penalty;
        end
        %% PENALTY CONTINUATION
        function obj = performPenaltyContinuation(obj)
            % Perform penalty continuation to increase the penalty factor
            % This method increases the penalty factor by a specified increment
            % until it reaches the maximum value defined in the penalty structure.
            % It checks if the current penalty factor is less than the maximum value,
            % and if so, it increases the penalty factor by the increment value.
            obj.m_penalty = min(obj.m_penalty + ...
                obj.m_penaltyStruct.inc,obj.m_penaltyStruct.max);
        end
        %% GET INTERPOLATION COEFFICIENT
        function interpCoeff = getInterpolationCoefficient(obj,rho)
            % Get the interpolation coefficient based on the density field
            % Inputs:
            %   rho: density field of size (numNodes x numScenarios)
            % Output:
            %   interpCoeff: interpolation coefficient of size (numNodes x numScenarios)    
            % This method computes the interpolation coefficient based on the density field.
            % It uses different formulas for different interpolation methods:
            % - For SIMP interpolation: E = ρᵖE₀, where ρ is the density and p is the penalty factor
            % - For RAMP interpolation: E =  ρ / (1 + q(1 − ρ)) E₀, where q is the penalty factor
            % - For linear interpolation: E = ρE₀
            % The method checks if the density field is specified and computes the interpolation coefficient accordingly.

            if (nargin < 2), error('Density field must be specified'); end
            if (any(rho(:) < 0) || any(rho(:) > 1))
                error('Density field must be between 0 and 1');
            end
            if strcmp(obj.m_interpolation,'simp')
                interpCoeff = max(obj.void_density_interpolated,(rho.^obj.m_penalty)); % for SIMP interpolation: E = ρᵖE₀

            elseif strcmp(obj.m_interpolation,'ramp')
                interpCoeff = max(1e-12,(rho./(1+obj.m_penalty*(1-rho)))) ; % for RAMP interpolation: E =  ρ / (1 + q(1 − ρ)) E₀, where q is the penalty factor

            else
                interpCoeff = max(obj.void_density,rho); % linear
            end
        end
        %% GET DERIVATIVE OF INTERPOLATION COEFFICIENT
        function interpCoeffGrad = getInterpolationCoefficientGrad(obj,rho)
            % Get the derivative of the interpolation coefficient
            % x is the density field, which should be between 0 and 1
            % Inputs:
            %   rho: density field (should be between 0 and 1)
            % Output:
            %   interpCoeffGrad: derivative of the interpolation coefficient
            % This method computes the derivative of the interpolation coefficient
            % based on the density field. It uses different formulas for different interpolation methods:
            % - For SIMP interpolation: E = ρᵖE₀, dE/dρ = p*ρ^(p-1)*E₀
            % - For RAMP interpolation: E =  ρ / (1 + q(1 − ρ)), dE/dρ = (1+q)/(1+q(1-ρ))^2 * E₀
            % - For linear interpolation: dE/dρ = 1
            % The method checks if the density field is specified and computes the derivative accordingly

            if (nargin < 2), error('Density field must be specified'); end

            if (any(rho(:) < 0) || any(rho(:) > 1))
                error('Density field must be between 0 and 1');
            end

            if strcmp(obj.m_interpolation,'simp')
                p = obj.m_penalty;
                interpCoeffGrad = p*rho.^(p-1); % for SIMP interpolation: E = ρᵖE₀, dE/dρ = p*ρ^(p-1)*E₀
            elseif strcmp(obj.m_interpolation,'ramp')
                q = obj.m_penalty;
                interpCoeffGrad = (1+q)./(1+q*(1-rho)).^2;  % for RAMP interpolation: E =  ρ / (1 + q(1 − ρ)) E₀, dE/dρ = (1+q)/(1+q(1-ρ))^2 * E₀
                interpCoeffGrad(rho<0.01) = 0.0;
            else
                interpCoeffGrad = 1; % linear
            end
        end
    end
end