%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discription:                                                              %
% This class implements a 2D fluid finite element analysis (FEA) solver     %
% using the finite element method (FEM) for incompressible flow problems.   %
% It is designed to solve the Navier-Stokes equations in a 2D domain with   %
% Dirichlet boundary conditions for velocity and Neumann boundary conditions%
% for pressure. The solver supports both steady and unsteady flow problems  %
% and includes features for handling Brinkman penalty factors, which are    %
% used to enforce boundary conditions and improve numerical stability.      %
% The class provides methods for assembling the global stiffness matrix,    %
% computing the residual vector, and solving the system of equations using  %
% a Newton-Raphson method. It also includes functionality for handling      %
% different boundary conditions, such as uniform and parabolic profiles,    %
% and for managing the degrees of freedom (DOFs) associated with velocity   %
% and pressure fields. The solver is designed to be flexible and can be     %
% adapted for various fluid flow problems by modifying the mesh, boundary.  %
% conditions, and material properties.                                      %
%                                                                           %
% This code is largely based on the MATLAB code written by  Joe             %
% Joe Alexandersen which is included in the:                                %
% 'utilities/thirdParty/topflow' directory.                                 %
% For more details, please refer to the original paper:                     %
% Alexandersen, Joe. "A detailed introduction to density-based              %
% topology optimisation of fluid flow problems with                         %
% implementation in MATLAB." Structural and Multidisciplinary               %
% Optimization 66.1 (2023): 12.                                             %
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

classdef fea2d_fluid < fea2d
    properties(GetAccess = 'public', SetAccess = 'private')

        m_fluxSegs;              % Segments for flux calculations
        m_internalForce = 0;     % Internal force variable, initialized to 0
        m_BCprofile;           % for each u or v segement, 0: uniform, 1: parabolic

        m_edofVecU;              % Vector of element degrees of freedom for velocity
        m_edofMatU;              % Matrix of element degrees of freedom for velocity
        m_edofVecP;              % Vector of element degrees of freedom for pressure
        m_edofMatP;              % Matrix of element degrees of freedom for pressure

        m_iJ;                    % Row indices for Jacobian matrix construction
        m_jJ;                    % Column indices for Jacobian matrix construction
        m_iR;                    % Row indices for residual vector construction
        m_jR;                    % Column indices for residual vector construction
        m_jE;
        m_freeDofsProjector;     % Projector matrix for free DOFs
        m_fixedDofsNullifier;    % Nullifier matrix for fixed DOFs
        m_residual;              % Residual vector

        m_velocity;              % Velocity field
        m_pressure;              % Pressure field

        m_qa0;                   % Initial value for Brinkman penalty factor
        m_qa;                    % Current value for Brinkman penalty factor
        m_alpha;                 % Current Brinkman penalty factor
        m_alphaGrad;             % Gradient of the Brinkman penalty factor
        m_alpha0;                % Initial Brinkman penalty factor
        m_alphaMin;              % Minimum Brinkman penalty factor
        m_alphaMax;              % Maximum Brinkman penalty factor
        m_qaStep;                % Current continuation step for Brinkman penalty factor

        % Newton solver parameters
        m_newtonTolerance = 1e-6; % Residual tolerance for Newton solver convergence
        m_maxNewtonIters = 25;    % Maximum number of Newton iterations

        m_continuationScheme;     % Heuristic continuation scheme
        m_maxNumContinuationSteps = 50; % Maximum number of continuation steps
        m_continuationSchemeSize; % Size of the continuation scheme
        m_numStreamlineSamples = 40;  % Number of samples for streamlines

        m_activeDesignDomain; % Active design domain, subset of existingElems (for example in wind tunnel)
        m_rectangleActiveDomainBbox; % design domain bbox from input
    end

    methods (Access = public)
        function obj = fea2d_fluid(brep,numElements,materials,...
                interpolation,numScenarios,penaltyStruct,uniformGrid)
            % set default values
            if nargin < 4, interpolation = 'none';end
            if nargin < 5, numScenarios = 1;end
            if nargin < 6; penaltyStruct = struct('min',1,'max',1,'inc',0);end
            if nargin < 7; uniformGrid = 0; end
            % construct
            numDOFperNode = 3; % [u,v,p]

            obj = obj@fea2d(brep,numElements,numDOFperNode,materials,...
                interpolation,numScenarios,penaltyStruct,uniformGrid); % call superclass

            obj.m_BCprofile = zeros(obj.m_numBndrySegs,obj.m_numDOFperNode,obj.m_numScenarios); % for each u or v segement, 0: uniform, 1: parabolic

            obj.m_velocity.u = zeros(obj.m_ny+1,obj.m_nx+1,obj.m_numScenarios);
            obj.m_velocity.v = zeros(obj.m_ny+1,obj.m_nx+1,obj.m_numScenarios);
            obj.m_velocity.norm = zeros(obj.m_ny+1,obj.m_nx+1,obj.m_numScenarios);
            obj.m_pressure = zeros(obj.m_ny+1,obj.m_nx+1,obj.m_numScenarios);

            % BRINKMAN PENALISATION
            obj.m_alphaMax= 2.5*obj.m_materials(1).mu/(0.01^2);
            obj.m_alphaMin = 2.5*obj.m_materials(1).mu/(100^2);
            obj.m_alpha0 = 2.5*obj.m_materials(1).mu/(0.1^2); % initial penalty value
            obj = obj.setupContinuationScheme();

            numNodes_y = obj.m_ny+1;
            numNodes_x = obj.m_nx+1;

            nodenrs = reshape(1:obj.m_numNodes,numNodes_y,numNodes_x);

            obj.m_edofVecU = reshape(2*nodenrs(1:end-1,1:end-1)+1,obj.m_numElems,1);
            obj.m_edofMatU = repmat(obj.m_edofVecU,1,8)+repmat([0 1 2*obj.m_ny+[2 3 0 1] -2 -1],obj.m_numElems,1);

            obj.m_edofVecP = reshape(nodenrs(1:end-1,1:end-1),obj.m_numElems,1);
            obj.m_edofMatP = repmat(obj.m_edofVecP,1,4)+repmat([1 obj.m_ny+[2 1] 0],obj.m_numElems,1);

            obj.m_edofMat = [obj.m_edofMatU 2*obj.m_numNodes+obj.m_edofMatP];
            obj.m_iJ = reshape(kron(obj.m_edofMat,ones(12,1))',144*obj.m_numElems,1);
            obj.m_jJ = reshape(kron(obj.m_edofMat,ones(1,12))',144*obj.m_numElems,1);
            obj.m_iR = reshape(obj.m_edofMat',12*obj.m_numElems,1);
            obj.m_jR = ones(12*obj.m_numElems,1);
            obj.m_jE = repmat(1:obj.m_numElems,12,1);

            obj.m_activeDesignDomain = obj.m_existingElems;
        end
        %% PENALTY CONTINUATION
        function obj = setAlphaValues(obj,alpha_min, alpha_max,alpha_0,qa)
            obj.m_alphaMax= alpha_max;
            obj.m_alphaMin = alpha_min;
            obj.m_alpha0 = alpha_0; % initial penalty value
            obj = obj.setupContinuationScheme(qa);
        end
        function obj = performPenaltyContinuation(obj)
            % Alexandersen, Joe. "A detailed introduction to density-based
            % topology optimisation of fluid flow problems with
            % implementation in MATLAB." Structural and Multidisciplinary
            % Optimization 66.1 (2023): 12.
            % Eq. 30 and 31
            if (obj.m_qaStep < obj.m_continuationSchemeSize)
                obj.m_qaStep = obj.m_qaStep + 1;
                obj.m_qa = obj.m_continuationScheme(qastep);
            end
        end
        function obj = set_qa(obj,qa)
            obj.m_qa = qa;
        end

        %% GET INTERPOLATION COEFFICIENT
        function obj = setupContinuationScheme(obj,qa)
            % Alexandersen, Joe. "A detailed introduction to density-based
            % topology optimisation of fluid flow problems with
            % implementation in MATLAB." Structural and Multidisciplinary
            % Optimization 66.1 (2023): 12.
            % Eq. 30 and 31
            idx = obj.m_activeDesignDomain==1;
            activeDesign = obj.m_design(idx);
            
            x0 = mean(activeDesign(:));
            obj.m_qa0 = (obj.m_alphaMax-obj.m_alpha0-x0*(obj.m_alphaMax-obj.m_alphaMin))/(x0*(obj.m_alpha0-obj.m_alphaMin));
            obj.m_continuationScheme = obj.m_qa0./[1 2 10 20]; % heurisitc
            obj.m_continuationSchemeSize = length(obj.m_continuationScheme);
            if nargin < 2
                obj.m_qa = obj.m_continuationScheme(1);
            else
                obj.m_qa = qa;
            end
        end

        function obj = computeInterpolationCoefficient(obj,xPhys)
            % interpolate Brinkman penalty factor
            % Alexandersen, Joe. "A detailed introduction to density-based
            % topology optimisation of fluid flow problems with
            % implementation in MATLAB." Structural and Multidisciplinary
            % Optimization 66.1 (2023): 12.
            % Eq. 8
            xPhys(obj.m_existingElems==0) = obj.void_density;
            obj.m_alpha = obj.m_alphaMin + ...
                (obj.m_alphaMax-obj.m_alphaMin)*(1-xPhys(:))./(1+obj.m_qa*xPhys(:));
        end
        %% GET DERIVATIVE OF INTERPOLATION COEFFICIENT
        function obj = computeInterpolationCoefficientGrad(obj,xPhys)
            alphaRange = obj.m_alphaMax - obj.m_alphaMin;
            dalpha_1 = (obj.m_qa*alphaRange*(xPhys(:) - 1))./(xPhys(:)*obj.m_qa + 1).^2;
            dalpha_2 = -alphaRange./(xPhys(:)*obj.m_qa + 1);
            obj.m_alphaGrad = dalpha_1 + dalpha_2;
        end
        %% MATERIAL PROPERTIES
        function obj = computeMaterialPropertiesMatrices(obj),return;end
        %% Internal Loads
        function obj = assembleInternalLoad(obj),return;end
        %% SOLVE
        function obj = solve(obj)
            % Code based on:
            % Alexandersen, Joe. "A detailed introduction to density-based
            % topology optimisation of fluid flow problems with
            % implementation in MATLAB." Structural and Multidisciplinary
            % Optimization 66.1 (2023): 12.
            obj.m_sol = zeros(obj.m_numDOFs,1);
            obj.m_sol(obj.m_fixedDOFs) = obj.m_fixed(obj.m_fixedDOFs);
            obj = obj.computeInterpolationCoefficient(obj.m_design);
            obj = obj.assembleK();
            %% NON-LINEAR NEWTON SOLVER
            iter = 0; fail = -1;
            while (fail ~= 1)
                iter = iter+1;
                % assemble residual
                obj = obj.assembleResidual(obj.m_sol);
                % evaluate residual norm
                if (iter == 1); r0 = norm(obj.m_residual); end
                r1 = norm(obj.m_residual);
                normR = r1/r0;
                if (normR < obj.m_newtonTolerance); break; end
                % assemble K
                obj = obj.assembleK();
                % calculate Newton step
                dS = -obj.m_K\obj.m_residual;
                % L2-norm line search
                Sp = obj.m_sol + 0.5*dS;
                obj = obj.assembleResidual(Sp);
                r2 = norm(obj.m_residual);
                Sp = obj.m_sol + 1.0*dS;
                obj = obj.assembleResidual(Sp);
                r3 = norm(obj.m_residual);
                % solution update with "optimal" damping
                lambda = max(0.01,min(1.0,(3*r1 + r3 - 4*r2)/(4*r1 + 4*r3 - 8*r2)));
                obj.m_sol = obj.m_sol + lambda*dS;
                % if fail, retry from zero solution
                if (iter == obj.m_maxNewtonIters && fail < 0)
                    iter = 0;
                    obj.m_sol(obj.m_freeDOFs) = 0.0;
                    normR=1; %#ok
                    fail = fail+1;
                end
                if (iter == obj.m_maxNewtonIters && fail < 1)
                    fail = fail+1;
                end
            end
        end
        function obj = assembleResidual(obj,S)
            dxv = obj.m_hx*ones(1,obj.m_numElems);
            dyv = obj.m_hy*ones(1,obj.m_numElems);
            muv = obj.m_materials(1).mu*ones(1,obj.m_numElems); % Assuming single fluid
            rhov = obj.m_materials(1).rho*ones(1,obj.m_numElems); % Assuming single fluid

            sR = RES(dxv,dyv,muv,rhov,obj.m_alpha(:)',S(obj.m_edofMat'));
            obj.m_residual = sparse(obj.m_iR,obj.m_jR,sR(:));
            obj.m_residual(obj.m_fixedDOFs) = 0;
        end
        %% ASSEMBLE STIFFNESS MATRIX
        function obj = assembleK(obj)
            dxv = obj.m_hx*ones(1,obj.m_numElems);
            dyv = obj.m_hy*ones(1,obj.m_numElems);
            muv = obj.m_materials(1).mu*ones(1,obj.m_numElems); % Assuming single fluid
            rhov = obj.m_materials(1).rho*ones(1,obj.m_numElems); % Assuming single fluid

            sJ = JAC(dxv,dyv,muv,rhov,obj.m_alpha(:)',obj.m_sol(obj.m_edofMat'));
            obj.m_K = sparse(obj.m_iJ,obj.m_jJ,sJ(:));
            obj.m_K = obj.m_fixedDofsNullifier'*obj.m_K*obj.m_fixedDofsNullifier ...
                + obj.m_freeDofsProjector;
        end
        %% BOUNDARY CONDITIONs
        function obj = fixUOfEdge(obj,boundaryEdges,uValue,profileMode,scenarioId)
            if (nargin < 4),profileMode=1;end %assuming parabolic
            if (nargin < 5),scenarioId=1;end
            obj.m_BCtype(boundaryEdges,1,scenarioId) = 1;
            obj.m_BCvalue(boundaryEdges,1,scenarioId) = uValue;
            obj.m_BCprofile(boundaryEdges,1,scenarioId) = profileMode;
        end
        function obj = fixVOfEdge(obj,boundaryEdges,vValue,profileMode,scenarioId)
            if (nargin < 4),profileMode=1;end %assuming parabolic
            if (nargin < 5),scenarioId=1;end
            obj.m_BCtype(boundaryEdges,2,scenarioId) = 1;
            obj.m_BCvalue(boundaryEdges,2,scenarioId) = vValue;
            obj.m_BCprofile(boundaryEdges,2,scenarioId) = profileMode;
        end
        function obj = fixPOfEdge(obj,boundaryEdges,pValue,scenarioId)
            if (nargin < 4),scenarioId=1;end
            obj.m_BCtype(boundaryEdges,3,scenarioId) = 1;
            obj.m_BCvalue(boundaryEdges,3,scenarioId) = pValue;
        end
        %% ASSEMBLE BOUNDARY CONDITIONS
        function obj = assembleBC(obj)
            nDOF = obj.m_numDOFs;
            obj.m_f = sparse(nDOF ,obj.m_numScenarios);
            for scenarioId = 1:obj.m_numScenarios
                isDirichlet = zeros(nDOF,1);
                dirValue = zeros(nDOF,1);
                isDirichlet(obj.m_fixedDOFs) = 1;
                for geomEdge = 1:size(obj.m_brep.segments,2)
                    typeu = obj.m_BCtype(geomEdge,1,scenarioId);
                    typev = obj.m_BCtype(geomEdge,2,scenarioId);
                    typep = obj.m_BCtype(geomEdge,3,scenarioId);
                    valueu = obj.m_BCvalue(geomEdge,1,scenarioId);
                    valuev = obj.m_BCvalue(geomEdge,2,scenarioId);
                    valuep = obj.m_BCvalue(geomEdge,3,scenarioId);
                    if (typeu == 1)
                        [nodes,points] = obj.findNodesOnEdge(geomEdge);

                        y_points = points(2,:);
                        y_length = max(y_points)-min(y_points);

                        if (y_length>1e-6) && (obj.m_BCprofile(geomEdge,1,scenarioId)==1)
                            y_normalized = (y_points-min(y_points))/y_length;
                            u = @(y) -4*y.^2+4*y;
                            U = valueu*u(y_normalized);
                        else
                            U = valueu;
                        end
                        obj.m_fixedNodes = [obj.m_fixedNodes;  nodes];
                        udof = 2*nodes-1;
                        isDirichlet(udof) = 1;
                        dirValue(udof) =  U;
                    end
                    if (typev == 1)
                        [nodes,points] = obj.findNodesOnEdge(geomEdge);

                        x_points = points(1,:);
                        x_length = max(x_points)-min(x_points);
                        if (x_length>1e-6) && (obj.m_BCprofile(geomEdge,2,scenarioId)==1)
                            x_normalized = (x_points-min(x_points))/x_length;
                            v = @(x) -4*x.^2+4*x;
                            V = valuev*v(x_normalized);
                        else
                            V = valuev;
                        end
                        obj.m_fixedNodes = [obj.m_fixedNodes;  nodes];
                        vdof = 2*nodes;
                        isDirichlet(vdof) = 1;
                        dirValue(vdof) =  V;
                    end
                    if (typep == 1)
                        nodes = obj.findNodesOnEdge(geomEdge);
                        obj.m_fixedNodes = [obj.m_fixedNodes;  nodes];
                        pdof = 2*(obj.m_nx+1)*(obj.m_ny+1) + nodes;
                        isDirichlet(pdof) = 1;
                        dirValue(pdof) =  valuep;
                    end
                end
                obj.m_fixed = dirValue;
                obj.m_fixedNodes = unique(obj.m_fixedNodes);
                dirichletDOF = find(isDirichlet(:) == 1);
                obj.m_fixedDOFs = unique([obj.m_fixedDOFs(:) dirichletDOF(:)]);
                allDOF = 1:obj.m_numDOFs;
                obj.m_freeDOFs  = setdiff(allDOF,obj.m_fixedDOFs);
            end
            % Nullspace matrices
            obj.m_freeDofsProjector = speye(obj.m_numDOFs); % Create an identity matrix of size doftot
            obj.m_fixedDofsNullifier = obj.m_freeDofsProjector ; % Copy FreeDofsProjector to FixedDofsNullifier
            obj.m_fixedDofsNullifier(obj.m_fixedDOFs, obj.m_fixedDOFs) = 0.0; % Set the entries of FixedDofsNullifier corresponding to fixed DOFs to 0
            obj.m_freeDofsProjector  = obj.m_freeDofsProjector - obj.m_fixedDofsNullifier; % Subtract FixedDofsNullifier from FreeDofsProjector
        end
        %% POST-PROCESSING
        function obj = postProcess(obj)
            % extract velocity and pressure fields
            % from the solution vector
            numNodes = (obj.m_nx+1)*(obj.m_ny+1);
            obj.m_velocity.u = zeros(obj.m_ny+1,obj.m_nx+1,obj.m_numScenarios);
            obj.m_velocity.v = zeros(obj.m_ny+1,obj.m_nx+1,obj.m_numScenarios);
            obj.m_velocity.norm = zeros(obj.m_ny+1,obj.m_nx+1,obj.m_numScenarios);
            obj.m_pressure = zeros(obj.m_ny+1,obj.m_nx+1,obj.m_numScenarios);

            for scenarioId = 1:obj.m_numScenarios
                obj.m_velocity.u(:,:,scenarioId) = reshape(obj.m_sol(1:2:2*numNodes),obj.m_ny+1,obj.m_nx+1);
                obj.m_velocity.v(:,:,scenarioId) = reshape(obj.m_sol(2:2:2*numNodes),obj.m_ny+1,obj.m_nx+1);
                obj.m_velocity.norm(:,:,scenarioId) = reshape(sqrt(obj.m_sol(1:2:2*numNodes).^2+obj.m_sol(2:2:2*numNodes).^2),obj.m_ny+1,obj.m_nx+1);
                obj.m_pressure(:,:,scenarioId) = reshape(obj.m_sol(2*numNodes+1:3*numNodes),obj.m_ny+1,obj.m_nx+1);
            end
        end

        function obj = createRectangularDesignDomain(obj,center,w,h)
            % setPseudoDensityInRectangle  Mark elements inside a rectangle with pseudo-density = 1
            %
            %   center : [cx, cy] center of rectangle (in global coordinates)
            %   w, h   : full width and height of rectangle
            %
            %   The rectangle is axis-aligned and elements are marked based on their
            %   center coordinates in obj.m_elemCoords. Only existing elements
            %   (m_existingElems ~= 0) are considered.

            cx = center(1);
            cy = center(2);

            xMinRect = cx - 0.5*w;
            xMaxRect = cx + 0.5*w;
            yMinRect = cy - 0.5*h;
            yMaxRect = cy + 0.5*h;

            obj.m_rectangleActiveDomainBbox = zeros(2,2);
            obj.m_rectangleActiveDomainBbox(1,1) = xMinRect;
            obj.m_rectangleActiveDomainBbox(1,2) = xMaxRect;
            obj.m_rectangleActiveDomainBbox(2,1) = yMinRect;
            obj.m_rectangleActiveDomainBbox(2,2) = yMaxRect;

            % Element center coordinates as column vectors
            xe = obj.m_elemCoords(1, :).';
            ye = obj.m_elemCoords(2, :).';

            % Mask of existing elements in linear indexing
            maskExisting = obj.m_existingElems(:) ~= 0;

            % Mask of elements whose centers are inside the rectangle
            insideRect = (xe >= xMinRect) & (xe <= xMaxRect) & ...
                (ye >= yMinRect) & (ye <= yMaxRect);

            % Combine masks and write into pseudo-density field
            idx = find(maskExisting & insideRect);
            obj.m_activeDesignDomain = zeros(size(obj.m_existingElems));
            obj.m_activeDesignDomain(idx) = 1;

            obj.m_design(idx) = 0.001;
        end
        %% RESULTS
        function obj = printFluidResults(obj)
            for scenarioId = 1:obj.m_numScenarios
                min_u = min(obj.m_velocity.u(:,:,scenarioId),[],'all');
                max_u = max(obj.m_velocity.u(:,:,scenarioId),[],'all');
                min_v = min(obj.m_velocity.v(:,:,scenarioId),[],'all');
                max_v = max(obj.m_velocity.v(:,:,scenarioId),[],'all');
                min_norm = min(obj.m_velocity.norm(:,:,scenarioId),[],'all');
                max_norm = max(obj.m_velocity.norm(:,:,scenarioId),[],'all');
                min_p = min(obj.m_pressure(:,:,scenarioId),[],'all');
                max_p = max(obj.m_pressure(:,:,scenarioId),[],'all');

                disp(['scenario: ' num2str(scenarioId) ' ' ...
                    'u: [' num2str(min_u) ' -> ' num2str(max_u) '], ' ...
                    'v: [' num2str(min_v) ' -> ' num2str(max_v) '], ' ...
                    'velocity norm: [' num2str(min_norm) ' -> ' num2str(max_norm) '], ' ...
                    'pressure: [' num2str(min_p) ' -> ' num2str(max_p) ']' ...
                    ]);
            end
        end
        %% PLOTTING
        function obj = plotVelocity(obj,method)
            % plot grid mesh, sets outside value to NaN, so they are
            % ignored in the plot
            plt = PlotId;
            cm = ColorMaps;
            if nargin < 2, method = 'SurfInterp'; end
            for scenarioId = 1:obj.m_numScenarios
                % Set up figure and initialize plot
                figure(plt.velocity + scenarioId); clf(gcf,'reset');
                set(gcf, 'Name', strjoin({'Velocity', num2str(scenarioId)}, ' '));
                % Reshape node coordinates for streamline and surface plots
                X = reshape(obj.m_nodeCoords(1,:), [obj.m_ny+1, obj.m_nx+1]);
                Y = reshape(obj.m_nodeCoords(2,:), [obj.m_ny+1, obj.m_nx+1]);
                F = obj.m_velocity.norm(:, :, scenarioId);
                F(obj.m_existingNodes == 0) = NaN;  % Mark non-existing nodes as NaN to hide them
                % Plot surface data using surf to get 3D context
                surf(X, Y, F, 'EdgeColor', 'none'); colormap(cm.velocity); hold on;
                view(2); colorbar;
                % Set up streamline starting points within the bounding box
                sx = obj.m_boundingBox(1, 1) * ones(1, obj.m_numStreamlineSamples);
                sy = linspace(obj.m_boundingBox(2, 1), obj.m_boundingBox(2, 2), obj.m_numStreamlineSamples);
                % Define the grid for the velocity data
                [XGrid, YGrid] = meshgrid(linspace(obj.m_boundingBox(1, 1), obj.m_boundingBox(1, 2), obj.m_nx + 1), ...
                    linspace(obj.m_boundingBox(2, 1), obj.m_boundingBox(2, 2), obj.m_ny + 1));
                % Reshape the velocity fields to fit the grid structure
                u = reshape(obj.m_velocity.u(:), obj.m_ny + 1, obj.m_nx + 1);
                v = reshape(obj.m_velocity.v(:), obj.m_ny + 1, obj.m_nx + 1);
                % Create streamlines using 'stream2' and extract streamline vertices
                streamData = stream2(XGrid, YGrid, u, -v, sx, sy);
                % Plot the streamlines manually, ensuring they are above the surface plot
                for k = 1:length(streamData)
                    streamlineCoords = streamData{k};
                    if ~isempty(streamlineCoords)
                        xData = streamlineCoords(:, 1);
                        yData = streamlineCoords(:, 2);
                        % Interpolate Z values from the surface F to project the streamline onto it
                        zData = interp2(X, Y, F, xData, yData, 'linear', NaN);  % Use NaN for out-of-bound points
                        % Plot 3D streamlines on top of the surface
                        plot3(xData, yData, zData, 'k', 'LineWidth', 2);  % Offset zData slightly to place it on top
                    end
                end

                % Set axis properties to match bounding box
                pbaspect(obj.m_boxSizes);
                axis on;
                xlim(obj.m_boundingBox(1, :));
                ylim(obj.m_boundingBox(2, :));

                if (strcmp(method,'VoxelModel') == 1)
                    grid on;
                elseif (strcmp(method,'SurfInterp') == 1)
                    hold on; grid off; shading interp;
                else
                    disp(['Method ' method 'for plotting is not implemented!']);
                end
            end
        end

        function obj = plotPressure(obj,method)
            % plot grid mesh, sets outside value to NaN, so they are
            % ignored in the plot
            plt = PlotId;
            cm = ColorMaps;
            if nargin < 2, method = 'SurfInterp'; end
            for scenarioId = 1:obj.m_numScenarios
                X = reshape(obj.m_nodeCoords(1,:),[obj.m_ny+1,obj.m_nx+1]);
                Y = reshape(obj.m_nodeCoords(2,:),[obj.m_ny+1,obj.m_nx+1]);
                F = obj.m_pressure(:,:,scenarioId);
                F(obj.m_existingNodes==0) = NaN;
                figure(plt.pressure+scenarioId); clf(gcf,'reset');
                set(gcf, 'Name', strjoin({'Pressure',num2str(scenarioId)},' ') );
                surf(X,Y,F); colormap(cm.pressure); view(2); colorbar;
                pbaspect(obj.m_boxSizes);axis on;
                xlim(obj.m_boundingBox(1,:));
                ylim(obj.m_boundingBox(2,:));
                if (strcmp(method,'VoxelModel') == 1)
                    grid on;
                elseif (strcmp(method,'SurfInterp') == 1)
                    grid off; shading interp;
                else
                    disp(['Method ' method 'for plotting is not implemented!']);
                end
            end
        end

        function obj = plotBoundaryCondition(obj)
            plt = PlotId;
            plbc = PlotBC;
            scale = 0.1*obj.m_modelScale;
            for scenarioId = 1:obj.m_numScenarios
                legend_exists = struct('noSlip_U', false, ...
                    'noSlip_V', false, ...
                    'flow_U', false, ...
                    'flow_V', false, ...
                    'fixed_P', false);
                legend_fields = {};
                legend_labels = {};

                fig = figure(plt.loading+scenarioId);clf(fig,'reset');
                % plot grid mesh, sets outside value to NaN, so they are ignored in the plot
                obj.plotGeometry(plt.loading+scenarioId,0);

                set(gcf, 'Name', strjoin({'Boundary Condition',num2str(scenarioId)},' ') );hold on;
                % mark all fixed nodes
                numNodes = (obj.m_nx+1)*(obj.m_ny+1);
                X = zeros(1,numNodes); Y = zeros(1,numNodes);
                hx = obj.m_boxSizes(1)/obj.m_nx;
                hy = obj.m_boxSizes(2)/obj.m_ny;
                for nodex = 1:obj.m_nx+1
                    for nodey = 1:obj.m_ny+1
                        nodeId = ((nodex-1)*(obj.m_ny+1) + nodey);
                        X(nodeId) = obj.m_boundingBox(1,1) + (nodex-1)*hx;
                        Y(nodeId) = obj.m_boundingBox(2,1) + (nodey-1)*hy;
                    end
                end

                for geomEdge = 1:size(obj.m_brep.segments,2)
                    typeu = obj.m_BCtype(geomEdge,1,scenarioId);
                    typev = obj.m_BCtype(geomEdge,2,scenarioId);
                    typep = obj.m_BCtype(geomEdge,3,scenarioId);
                    valueu = obj.m_BCvalue(geomEdge,1,scenarioId);
                    valuev = obj.m_BCvalue(geomEdge,2,scenarioId);
                    [~,points] = obj.findNodesOnEdge(geomEdge);
                    if isempty(points),continue;end
                    x_points = points(1,:);
                    y_points = points(2,:);
                    if (typeu == 1) % flow
                        y_length = max(y_points)-min(y_points);
                        if y_length>1e-6 && valueu > 1e-6

                            y_normalized = (y_points-min(y_points))/y_length;
                            u = @(y) -4*y.^2+4*y;
                            U = u(y_normalized);
                            if (obj.m_BCprofile(geomEdge,1,scenarioId)==0)
                                U = ones(size(U));
                            end
                            start = [x_points; y_points];
                            stop = start;
                            stop(1,:) = start(1,:) + scale*U;
                            obj.drawArrow(start',stop',plbc.flow_U.color,2);hold on;
                            % legend
                            if ~legend_exists.flow_U
                                flow_U = plot(NaN, NaN, plbc.flow_U.marker, ...
                                    'MarkerEdgeColor', plbc.flow_U.color, ...
                                    'MarkerFaceColor', plbc.flow_U.color); hold on;
                                legend_fields = [legend_fields;flow_U ]; %#ok
                                legend_labels = [legend_labels,'flow $u$']; %#ok
                                legend_exists.flow_U = true;
                            end
                        else % no-slip
                            scatter(x_points,y_points,'MarkerEdgeColor',plbc.noSlip_U.color,'Marker',plbc.noSlip_U.marker);hold on;
                            % legend
                            if ~legend_exists.noSlip_U
                                noSlip_U = plot(NaN, NaN, ...
                                    plbc.noSlip_U.marker, ...
                                    'MarkerEdgeColor', plbc.noSlip_U.color); hold on;
                                legend_fields = [legend_fields;noSlip_U ]; %#ok
                                legend_labels = [legend_labels,'no-slip $u$']; %#ok
                                legend_exists.noSlip_U = true;
                            end
                        end
                    end

                    if (typev == 1)
                        x_length = max(x_points)-min(x_points);
                        if x_length>1e-6 && valuev > 1e-6
                            x_normalized = (x_points-min(x_points))/x_length;
                            v = @(x) -4*x.^2+4*x;
                            V = v(x_normalized);
                            if (obj.m_BCprofile(geomEdge,2,scenarioId)==0)
                                V = ones(size(V));
                            end
                            start = [x_points; y_points];
                            stop = start;
                            stop(2,:) = start(2,:) + scale*V;
                            obj.drawArrow(start',stop',plbc.flow_V.color,2);hold on;
                            % legend
                            if ~legend_exists.flow_V
                                flow_V = plot(NaN, NaN,plbc.flow_V.marker,'MarkerEdgeColor',plbc.flow_V.color,'MarkerFaceColor', plbc.flow_V.color);
                                legend_fields = [legend_fields;flow_V ]; %#ok
                                legend_labels = [legend_labels,'flow $v$']; %#ok
                                legend_exists.flow_V = true;
                            end
                        else
                            scatter(x_points,y_points,'MarkerEdgeColor',plbc.noSlip_V.color,'Marker',plbc.noSlip_V.marker);hold on;
                            % legend
                            if ~legend_exists.noSlip_V
                                noSlip_V = plot(NaN, NaN, plbc.noSlip_V.marker, 'MarkerEdgeColor', plbc.noSlip_V.color);
                                legend_fields = [legend_fields;noSlip_V ]; %#ok
                                legend_labels = [legend_labels,'no-slip $v$']; %#ok
                                legend_exists.noSlip_V = true;
                            end
                        end

                    end

                    if (typep == 1)
                        scatter(x_points,y_points,'Marker',plbc.fixed_P.marker,'MarkerFaceColor',plbc.fixed_P.color);hold on;
                        % legend
                        if ~legend_exists.fixed_P
                            fixed_P = plot(NaN, NaN,plbc.fixed_P.marker,'MarkerFaceColor',plbc.fixed_P.color);
                            legend_fields = [legend_fields;fixed_P ]; %#ok
                            legend_labels = [legend_labels,'fixed $p$']; %#ok
                            legend_exists.fixed_P = true;
                        end
                    end
                end


                % plot active domain (if is specified)
                if ~isempty(obj.m_rectangleActiveDomainBbox)
                    x = obj.m_rectangleActiveDomainBbox(1,[1 2 2 1]);
                    y = obj.m_rectangleActiveDomainBbox(2,[1 1 2 2]);

                    hDom = patch(x, y, [0.8 0.8 0.8], ...
                        'EdgeColor','k', ...
                        'FaceAlpha',0.4, ...
                        'DisplayName','Design domain');
                    legend_fields = [legend_fields; hDom(:)]; %#ok
                    legend_labels = [legend_labels, {'active design domain'}]; %#ok
                end


                legend(legend_fields,legend_labels, ...
                    'Location', 'northeastoutside');

                pbaspect(obj.m_boxSizes);axis on;axis tight;
            end
        end
    end
end