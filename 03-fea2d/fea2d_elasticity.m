%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                              %
% This class implements the finite element analysis for 2D elasticity       %
% problems. It inherits from the fea2d class and provides methods to        %
% compute the elasticity matrix, stiffness matrix, and perform finite       %
% element analysis for elasticity problems.                                 %
%                                                                           %
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

classdef fea2d_elasticity < fea2d
    properties(GetAccess = 'public', SetAccess = 'protected')
        m_acceleration; % acceleration
        m_D; % elasticity matrix
        m_def; % deformation vector
        m_strainTensor;  % vector of strain tensors at each element     
        m_stressTensor;   % vector of stress tensors at each element
        m_vonMisesElems; % vector of von Mises stresses at each element
        m_vonMisesNodes;  % vector of von Mises stresses at each node
        m_principalStress;  % principal stresses at each element along x and y axes, as well as tension and compression
        m_maxDef;
        m_maxStress;
    end
    methods (Access = public)
        %% CONSTRUCT
        function obj = fea2d_elasticity(brep,numElements,materials,...
                vectorize,numScenarios,interpolation,penaltyStruct,uniformGrid)
            % set default values
            if nargin < 4; vectorize = true;end
            if nargin < 5, numScenarios = 1;end
            if nargin < 6, interpolation = 'none';end
            if nargin < 7; penaltyStruct = struct('min',1,'max',1,'inc',0);end
            if nargin < 8; uniformGrid = 0;end
            % construct
            numDOFperNode = 2; % u v
            obj = obj@fea2d(brep,numElements,numDOFperNode,materials,...
                interpolation,numScenarios,penaltyStruct,uniformGrid); % call superclass
            obj.m_materialIndices = ones(obj.m_ny,obj.m_nx);
            obj.m_acceleration = zeros(numScenarios,2);
            obj.m_principalStress = struct('x',[],'y',[], ...
                'tension',[],'compression',[]);
            %%
            obj.m_vectorize = vectorize;
            if obj.m_vectorize 
                %% PREPARE FINITE ELEMENT ANALYSIS
                nelx = obj.m_nx;
                nely = obj.m_ny;
                nodenrs = reshape(1:(nelx+1)*(nely+1),nely+1,nelx+1);
                obj.m_edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
                obj.m_edofMat = repmat(obj.m_edofVec,1,8)+ ...
                    repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
                edofMat2dofOrder = [7 5 3 1 8 6 4 2];
                obj.m_edofMat = obj.m_edofMat(:,edofMat2dofOrder);
                % The Kronecker product (kron) is a block matrix where
                % each element of the first matrix is multiplied by the
                % entire second matrix.
                obj.m_iK = reshape(kron(obj.m_edofMat,ones(8,1))',64*nelx*nely,1);
                obj.m_jK = reshape(kron(obj.m_edofMat,ones(1,8))',64*nelx*nely,1);
            end
        end
        %% COMPUTE ELASTICITY MATRIX
        function obj = computeMaterialPropertiesMatrices(obj)
            for matId = 1:obj.m_numMaterials
                % isotropic
                E = obj.m_materials(matId).E;
                nu = obj.m_materials(matId).nu;
                obj.m_D{matId} = E/(1-nu^2)*[1 nu 0; nu 1 0;0 0 (1-nu)/2];
            end
            % compute template stiffness matrix
            obj.m_KE = obj.integrateKOverElemQ4();
        end
        function D = getElasticityMatrix(obj,matId)
            % get elasticity matrix for a given material
            if nargin == 1
                D = obj.m_D;
            else
                D = obj.m_D{matId};
            end
        end
        %% ASSEMBLE STIFFNESS MATRIX
        function obj = assembleK(obj)
            % Assemble the global stiffness matrix
            % The global stiffness matrix is assembled by integrating the
            % elasticity matrix over each element and summing the contributions
            % from all elements.
            % The global stiffness matrix is a sparse matrix, and the
            % contributions from each element are stored in triplet form
            % (RowTriplets, ColTriplets, EntryTriplets).
            % The global stiffness matrix is then constructed using the
            % sparse function, which creates a sparse matrix from the triplet
            % form.
            % The global stiffness matrix is symmetric, so it is averaged with
            % its transpose to ensure symmetry.
            % The global stiffness matrix is stored in the m_K property of the
            % fea2d_elasticity object.
            % The global stiffness matrix is used in the finite element analysis
            % to solve for the displacements and stresses in the structure.
            nelx = obj.m_nx;
            nely = obj.m_ny;
            if ~obj.m_vectorize 
                nDOF = 2*(nelx+1)*(nely+1);
                nElements = nelx*nely;
                nzmax = 64*nElements;
                RowTriplets = zeros(nzmax,1);
                ColTriplets = zeros(nzmax,1);
                EntryTriplets = zeros(nzmax,1);
                index = 1;
                for elx = 1:nelx
                    for ely = 1:nely
                        nodes = [((elx-1)*(nely+1) + ely),...
                            ((elx)*(nely+1) + ely),...
                            ((elx)*(nely+1) + ely+1),...
                            ((elx-1)*(nely+1) + ely+1)];
                        rho = obj.m_design(ely,elx);
                        interpCoeff = obj.getInterpolationCoefficient(rho);     
                        KElem = interpCoeff * obj.m_KE;     
                        dof = [2*nodes-1; 2*nodes];
                        dof = reshape(dof,1,8);
                        temp = dof(ones(1,8),:);
                        colIndex = reshape(temp',1,8^2);    
                        rowIndex = reshape(temp,1,8^2);
                        entries = reshape(KElem',1,8^2);
                        RowTriplets(index:index+8^2-1,1) = rowIndex';
                        ColTriplets(index:index+8^2-1,1) = colIndex';
                        EntryTriplets(index:index+8^2-1,1) = entries';
                        index = index+8^2;
                    end
                end
                obj.m_K = sparse(RowTriplets,ColTriplets,EntryTriplets,nDOF,nDOF);
            else
                nelx = obj.m_nx;        nely = obj.m_ny;
                xPhys = obj.getInterpolationCoefficient(obj.m_design);
                sK = reshape(obj.m_KE(:)*xPhys(:)',64*nelx*nely,1);
                obj.m_K = sparse(obj.m_iK,obj.m_jK,sK);
                obj.m_K  = (obj.m_K +obj.m_K')/2;
            end
        end
        %% ASSEMBLE BOUNDARY CONDITIONS
        function obj = assembleInternalLoad(obj)
            % Assemble the internal loads for the finite element analysis.
            % The internal loads are computed by integrating the body force
            % over each element and summing the contributions from all elements.
            % The internal loads are stored in the m_fBody property of the
            % fea2d_elasticity object.
            obj.m_fE = obj.integrateBodyForceOverElemQ4();
            obj.m_fBody = zeros(obj.m_numDOFs,obj.m_numScenarios);
            nelx = obj.m_nx;
            nely = obj.m_ny;
            %% Assemble internal loads
            for scenarioId = 1:obj.m_numScenarios
                for elx = 1:nelx
                    for ely = 1:nely
                        nodes = [((elx-1)*(nely+1) + ely),...
                            ((elx)*(nely+1) + ely),...
                            ((elx)*(nely+1) + ely+1),...
                            ((elx-1)*(nely+1) + ely+1)];
                        dof = [2*nodes-1; 2*nodes];
                        dof = reshape(dof,8,1);
                        fbody = obj.m_fE(scenarioId,:);

                        obj.m_fBody(dof,scenarioId) = obj.m_fBody(dof,scenarioId) + ...
                            obj.m_design(ely,elx)*fbody';
                    end
                end
            end
        end

        function obj = assembleBC(obj)
            % Assemble the boundary conditions for the finite element analysis.
            % The boundary conditions are computed by integrating the surface
            % forces (Neumann data) over the boundary edges and applying the
            % Dirichlet conditions (fixed displacements) to the nodes.
            % The boundary conditions are stored in the m_f property of the
            % fea2d_elasticity object.
            % The m_fixedDOFs property contains the degrees of freedom that are
            % fixed (Dirichlet conditions), and the m_freeDOFs property contains
            % the degrees of freedom that are free (not fixed).
            % The m_forcedNodes property contains the nodes that are subjected to
            % surface forces (Neumann data).
            % The m_fixed property contains the values of the fixed displacements
            % at the fixed nodes.

            [xi_GQ, wt_GQ] = obj.GaussQLine();
            N1D = cell(1,length(xi_GQ));
            for i = 1:length(xi_GQ)
                N1D{i} = obj.edgeShapeFunction(xi_GQ(i));
            end
            nDOF = obj.m_numDOFs;
            obj.m_f = sparse(nDOF ,obj.m_numScenarios);

            obj.m_forcedNodes = [];
            obj.m_forcedNodes{obj.m_numScenarios} = [];
            for scenarioId = 1:obj.m_numScenarios
                % Assemble surface force (Neumann data)
                fBoundary = zeros(nDOF,1);
                isDirichlet = zeros(nDOF,1);
                dirValue = zeros(nDOF,1);
                isDirichlet(obj.m_fixedDOFs) = 1;

                for geomEdge = 1:size(obj.m_brep.segments,2)
                    typeu = obj.m_BCtype(geomEdge,1,scenarioId);
                    typev = obj.m_BCtype(geomEdge,2,scenarioId);
                    valueu = obj.m_BCvalue(geomEdge,1,scenarioId);
                    valuev = obj.m_BCvalue(geomEdge,2,scenarioId);
                    if (typeu == 0 ) && (abs(valueu) > 0)% force along x
                        boundarySegments = find(obj.m_edges(5,:) == geomEdge);
                        for seg = boundarySegments
                            nodes = obj.m_edges(1:obj.m_nodesPerEdge,seg);
                            obj.m_forcedNodes{scenarioId} = unique([obj.m_forcedNodes{scenarioId} ;nodes]);
                            udof = 2*nodes-1;
                            fBoundaryElem = obj.integrateOverBoundary(geomEdge,seg,wt_GQ,N1D,1,scenarioId);
                            fBoundary(udof) = fBoundary(udof) + fBoundaryElem;
                        end
                    end
                    if (typev == 0 ) && (abs(valuev) > 0)% force along y

                        boundarySegments = find(obj.m_edges(5,:) == geomEdge);
                        for seg = boundarySegments
                            nodes = obj.m_edges(1:obj.m_nodesPerEdge,seg);
                            obj.m_forcedNodes{scenarioId} = unique([obj.m_forcedNodes{scenarioId} ;nodes]);
                            vdof = 2*nodes;
                            fBoundaryElem = obj.integrateOverBoundary(geomEdge,seg,wt_GQ,N1D,2,scenarioId);

                            fBoundary(vdof) = fBoundary(vdof) + fBoundaryElem;
                        end
                    end
                    if (typeu == 1) % Dirichlet in x
                        boundarySegments = find(obj.m_edges(5,:) == geomEdge);
                        for seg = boundarySegments
                            nodes = obj.m_edges(1:obj.m_nodesPerEdge,seg);
                            nodes = unique(nodes(:));
                            obj.m_fixedNodes = [obj.m_fixedNodes;  nodes];
                            udof = 2*nodes-1;
                            isDirichlet(udof) = 1;
                            dirValue(udof) =  valueu;
                        end
                    end
                    if (typev == 1)  % Dirichlet in y
                        boundarySegments = find(obj.m_edges(5,:) == geomEdge);
                        for seg = boundarySegments
                            nodes = obj.m_edges(1:obj.m_nodesPerEdge,seg);
                            nodes = unique(nodes(:));
                            obj.m_fixedNodes = [obj.m_fixedNodes;  nodes];
                            vdof = 2*nodes;
                            isDirichlet(vdof) = 1;
                            dirValue(vdof) =  valuev;
                        end
                    end
                end
                obj.m_fixed = dirValue;
                obj.m_fixedNodes = unique(obj.m_fixedNodes);
                dirichletDOF = find(isDirichlet(:) == 1);
                obj.m_fixedDOFs = unique([obj.m_fixedDOFs(:) dirichletDOF(:)]);
                allDOF = 1:obj.m_numDOFs;
                obj.m_freeDOFs  = setdiff(allDOF,obj.m_fixedDOFs);

                obj.m_f(:,scenarioId) = obj.m_f(:,scenarioId) + fBoundary;
            end
        end


        function fBoundaryElem = integrateOverBoundary(obj,geomEdge,seg,wt_GQ,N1D,dof,scenarioId)
            % Integrate the boundary condition over the edge segment
            % Inputs:
            %   geomEdge: the geometric edge index
            %   seg: the segment index within the geometric edge            
            %   wt_GQ: the Gauss quadrature weights
            %   N1D: the shape functions for the 1D element
            %   dof: the degree of freedom (1 for u, 2 for v)    
            %   scenarioId: the scenario index            
            % Outputs:
            %   fBoundaryElem: the force vector for the boundary element            
            % This method integrates the boundary condition over the edge segment and returns the force vector for the boundary element.

            nodes = obj.m_edges(1:2,seg);
            xNodes = obj.m_nodeCoords(1,nodes);
            yNodes = obj.m_nodeCoords(2,nodes);
            dx = xNodes(2)-xNodes(1);
            dy = yNodes(2)-yNodes(1);
            L = sqrt(dx^2 + dy^2);
            vec = [dx dy 0]/L;
            zVec = [0 0 1];
            normal = cross(vec,zVec);
            nx = normal(1); %#ok<NASGU>
            ny = normal(2); %#ok<NASGU>
            fBoundaryElem = zeros(numel(N1D{1}),1);
            for g = 1:length(wt_GQ)
                N = N1D{g};
                f =  obj.m_BCvalue(geomEdge,dof,scenarioId); % dof is either 1 (u) or 2 (v)
                fBoundaryElem = fBoundaryElem + wt_GQ(g)*(L/2)*N*f;
            end
        end
        %% POST-PROCESSING
        function obj = postProcess(obj,principalStressesFlag)
            
            if nargin == 1, principalStressesFlag = false;end
            obj = obj.computeDeformation();
            obj = obj.computeElemStresses();

            obj.m_vonMisesNodes = zeros(obj.m_ny+1,obj.m_nx+1,obj.m_numScenarios);
            obj.m_maxDef = zeros(1,obj.m_numScenarios);
            obj.m_maxStress = zeros(1,obj.m_numScenarios);
            for scenarioId = 1:obj.m_numScenarios
                [obj,obj.m_vonMisesNodes(:,:,scenarioId)] = obj.computeNodalField(obj.m_vonMisesElems(:,:,scenarioId));

                obj.m_maxDef(scenarioId) = max(obj.m_def(:,:,scenarioId),[],'all');
                obj.m_maxStress(scenarioId) = max(obj.m_vonMisesElems(:,:,scenarioId),[],'all');
            end
            if principalStressesFlag,obj = obj.computePrincipalStress();end
        end
        %% BOUNDARY/LOADING CONDITIONS
        function obj = fixXOfEdge(obj,boundaryEdges)
            obj.m_BCtype(boundaryEdges,1,1) = 1;
        end
        function obj = fixYOfEdge(obj,boundaryEdges)
            obj.m_BCtype(boundaryEdges,2,1) = 1;
        end
        function obj = applyXForceOnEdge(obj,boundaryEdges,force,scenarioId)
            if (nargin < 4),scenarioId=1;end
            % convert force into pressure
            obj.m_BCvalue(boundaryEdges,1,scenarioId) = force./obj.m_segLengths(boundaryEdges);
        end
        function obj = applyYForceOnEdge(obj,boundaryEdges,force,scenarioId)
            if (nargin < 4),scenarioId=1;end
            % convert force into pressure
            obj.m_BCvalue(boundaryEdges,2,scenarioId) = force./obj.m_segLengths(boundaryEdges);
        end
        function obj = applyAcceleration(obj,value,scenarioId)
            if nargin < 3, scenarioId=1;end
            obj.m_acceleration(scenarioId,:) = value;
        end

        %% RESULTS
        function obj = printElascticityResults(obj)
            for scenarioId = 1:obj.m_numScenarios
                disp([ 'scenario: ' num2str(scenarioId) ...
                    ', max. Deformation: ' num2str(obj.m_maxDef(scenarioId)) ...
                    ', max. vonMises: ' num2str(obj.m_maxStress(scenarioId))]);
            end
        end
        %% PLOTTING
        function obj= plotBoundaryCondition(obj)
            plt = PlotId;
            plbc = PlotBC;
            for scenarioId = 1:obj.m_numScenarios
                legend_fields = {};
                legend_labels = {};
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

                index = (rem(obj.m_fixedDOFs,2)==1);
                fixedXNodes = (obj.m_fixedDOFs(index)-1)/2+1;
                X_xFixed = X(fixedXNodes);
                Y_xFixed = Y(fixedXNodes);
                if ~isempty(fixedXNodes)
                    fixed_U = plot(X_xFixed',Y_xFixed',plbc.fixed_U.marker,'MarkerEdgeColor', plbc.fixed_U.color);
                    hold on;
                    % legend
                    legend_fields = [legend_fields;fixed_U ]; %#ok
                    legend_labels = [legend_labels,'fixed $u$']; %#ok
                end
                index = (rem(obj.m_fixedDOFs,2)==0);
                fixedYNodes = (obj.m_fixedDOFs(index))/2;
                X_yFixed = X(fixedYNodes);
                Y_yFixed = Y(fixedYNodes);
                if ~isempty(fixedYNodes)
                    fixed_V = plot(X_yFixed',Y_yFixed',plbc.fixed_V.marker,'MarkerEdgeColor', plbc.fixed_V.color);
                    hold on;
                    % legend
                    legend_fields = [legend_fields;fixed_V ]; %#ok
                    legend_labels = [legend_labels,'fixed $v$']; %#ok
                end
                % mark all forced nodes
                nodes = obj.m_forcedNodes{scenarioId};
                if ~isempty(nodes)
                    % nodes = nodes(1:2:end);
                    X_forced = X(nodes);
                    Y_forced = Y(nodes);
                    Fx = obj.m_f(2*nodes-1,scenarioId);
                    Fy = obj.m_f(2*nodes,scenarioId);
                    scale = 0.1*obj.m_modelScale;
                    normF = sqrt(Fx.^2 + Fy.^2);
                    if (normF > 0)
                        Fx = Fx./normF;
                        Fy = Fy./normF;
                        Fx = Fx(:)';
                        Fy = Fy(:)';
                        plot(X_forced, Y_forced, 'd', 'MarkerFaceColor', plbc.force.color, 'LineStyle', 'none'); hold on;
                        start = [X_forced; Y_forced];
                        stop = start + scale*[Fx; Fy];
                        obj.drawArrow(start',stop',plbc.force.color,2);hold on;
                        % legend
                        force = plot(NaN, NaN,plbc.force.marker,'MarkerEdgeColor',plbc.force.color,'MarkerFaceColor', plbc.force.color); %hold on
                        legend_fields = [legend_fields;force ]; %#ok
                        legend_labels = [legend_labels,'force']; %#ok
                    end
                end
                % acceleration
                bx = obj.m_acceleration(scenarioId,1);
                by = obj.m_acceleration(scenarioId,2);
                bNorm = sqrt(bx^2+by^2);
                if bNorm > 1e-6
                    % Get axis limits
                    x_limits = xlim;
                    y_limits = ylim;
                    scale = 0.15*obj.m_modelScale;
                    start = obj.findPointinEmptyRegion(x_limits,y_limits,[[X_xFixed,X_yFixed];[Y_xFixed,Y_yFixed]]);
                    stop = start + scale*[bx/bNorm; by/bNorm];
                    obj.drawArrow(start',stop',plbc.acceleration.color,4);hold on;
                    % legend
                    acceleration = plot(NaN, NaN,plbc.acceleration.marker, ...
                        'MarkerEdgeColor',plbc.acceleration.color, ...
                        'MarkerFaceColor', plbc.acceleration.color); hold on
                    legend_fields = [legend_fields;acceleration ]; %#ok
                    legend_labels = [legend_labels,'body force']; %#ok
                end
                legend(legend_fields,legend_labels, ...
                    'Location', 'northeastoutside');
                pbaspect(obj.m_boxSizes);axis on;axis tight;
            end
        end
        function obj = plotDeformation(obj,shadingType)
            % plot grid mesh, sets outside value to NaN, so they are ignored in the plot
            plt = PlotId;
            cm = ColorMaps;
            if nargin == 1, shadingType = 'interp'; end
            scale = 0.05*obj.m_modelScale/max(obj.m_def(:));
            for scenarioId = 1:obj.m_numScenarios
                X = reshape(obj.m_nodeCoords(1,:),[obj.m_ny+1,obj.m_nx+1]);
                Y = reshape(obj.m_nodeCoords(2,:),[obj.m_ny+1,obj.m_nx+1]);

                delX = reshape(obj.m_sol(1:2:end,scenarioId),[obj.m_ny+1,obj.m_nx+1]);
                delY = reshape(obj.m_sol(2:2:end,scenarioId),[obj.m_ny+1,obj.m_nx+1]);

                X = X + scale*delX; % add deformation to mesh
                Y = Y + scale*delY;
                F = obj.m_def(:,:,scenarioId); F(obj.m_solidNodes==0) = NaN;
                figure(plt.deformation+scenarioId);
                set(gcf, 'Name', strjoin({'Deformation',num2str(scenarioId)},' ') );
                pdegplot(obj.m_pdeGeom);hold on
                surf(X,Y,F); colormap(cm.deformation); view(2); colorbar;
                pbaspect(obj.m_boxSizes);axis off;
                shading(shadingType);

                if (strcmp(shadingType,'faceted') == 1)
                    obj.plotWireMesh(plt.deformation+scenarioId);
                end
            end
        end
        function obj = plotVonMisesStress(obj,shadingType)
            % plot grid mesh, sets outside value to NaN, so they are
            % ignored in the plot
            plt = PlotId;
            cm = ColorMaps;
            if nargin == 1, shadingType = 'interp'; end
            for scenarioId = 1:obj.m_numScenarios
                X = reshape(obj.m_nodeCoords(1,:),[obj.m_ny+1,obj.m_nx+1]);
                Y = reshape(obj.m_nodeCoords(2,:),[obj.m_ny+1,obj.m_nx+1]);
                F = obj.m_vonMisesNodes(:,:,scenarioId);
                Fe = obj.m_vonMisesElems(:,:,scenarioId);

                % Nodal vs. elemental stress can vary a lot in range
                % Since all plots are nodal, we uniformly scale the plot
                % so the colorbar values are consistent
                maxFn = max(F(:));
                if (maxFn > 1e-6)
                    scale = max(Fe(:))/maxFn;
                    F = scale*F;
                end

                F(obj.m_solidNodes==0) = NaN;
                figure(plt.von_mises+scenarioId);
                set(gcf, 'Name', strjoin({'von Mises Stress ',num2str(scenarioId)},' ') );
                surf(X,Y,F); colormap(cm.von_mises); view(2);colorbar;
                pbaspect(obj.m_boxSizes);axis off;
                xlim(obj.m_boundingBox(1,:));
                ylim(obj.m_boundingBox(2,:));
                shading(shadingType);
            end
        end
        function obj = plotPrincipalStress(obj,method)
            compactness = 2;
            plt = PlotId;
            if (isempty(obj.m_principalStress.x)), obj = obj.computePrincipalStress(); end
            if nargin == 1, method = 'StreamLine'; end
            for scenarioId = 1:obj.m_numScenarios
                figure(plt.principal_stress+scenarioId);
                pdegplot(obj.m_pdeGeom);hold on
                set(gcf, 'Name', strjoin({'Principal Stress',num2str(scenarioId)},' ') );
                if (strcmp(method,'StreamLine') == 1)
                    X = reshape(obj.m_elemCoords(1,:),[obj.m_ny,obj.m_nx]);
                    Y = reshape(obj.m_elemCoords(2,:),[obj.m_ny,obj.m_nx]);
                    even_stream_line(X, Y, ...
                        obj.m_principalStress.tension.u(:,:,scenarioId).*(obj.m_solidElems==1), ...
                        obj.m_principalStress.tension.v(:,:,scenarioId).*(obj.m_solidElems==1), compactness,4, 'Color', 'r'); hold on

                    even_stream_line(X, Y, ...
                        obj.m_principalStress.compression.u(:,:,scenarioId).*(obj.m_solidElems==1), ...
                        obj.m_principalStress.compression.v(:,:,scenarioId).*(obj.m_solidElems==1), compactness,4,  'Color', 'b');

                else
                    quiver(obj.m_principalStress.x(:,:,scenarioId), ...
                        obj.m_principalStress.y(:,:,scenarioId), ...
                        obj.m_principalStress.tension.u(:,:,scenarioId).*(obj.m_solidElems==1), ...
                        obj.m_principalStress.tension.v(:,:,scenarioId).*(obj.m_solidElems==1),'ShowArrowHead','off','Color','r'); hold on

                    quiver(obj.m_principalStress.x(:,:,scenarioId), ...
                        obj.m_principalStress.y(:,:,scenarioId), ...
                        obj.m_principalStress.compression.u(:,:,scenarioId).*(obj.m_solidElems==1), ...
                        obj.m_principalStress.compression.v(:,:,scenarioId).*(obj.m_solidElems==1),'ShowArrowHead','off','Color','b');

                end
                pbaspect(obj.m_boxSizes);axis off; axis tight
                xlim(obj.m_boundingBox(1,:));
                ylim(obj.m_boundingBox(2,:));
                legend('Tension', 'Compression', ...
                    'Location', 'northeastoutside');
            end
        end
        %% INTEGRATE STIFFNESS MATRIX AT QUADRILATERAL ELEMENT
        function KE = integrateKOverElemQ4(obj)
            [xi_GQ,eta_GQ,wt_GQ]= obj.GaussQuad();
            NCell = cell(1,length(xi_GQ));
            gradNCell = cell(1,length(xi_GQ));
            for i = 1:length(xi_GQ)
                [NCell{i},gradNCell{i}] = obj.QuadShapeFunction(xi_GQ(i),eta_GQ(i));
            end
            KE = zeros(8,8,obj.m_numMaterials);
            Z = zeros(1,4);
            xNodes = obj.m_hx*[0,1,1,0];
            yNodes = obj.m_hy*[0,0,1,1];
            for matId = 1:obj.m_numMaterials
                for g = 1:length(wt_GQ)
                    gradN = gradNCell{g};
                    J = obj.Jacobian(xNodes,yNodes,xi_GQ(g),eta_GQ(g));
                    dJ = det(J);
                    T = J'\gradN;
                    B = [T(1,:) Z; Z T(2,:); T(2,:) T(1,:)];
                    KE(:,:,matId) = KE(:,:,matId) + wt_GQ(g)*dJ*B'*obj.m_D{matId}*B;
                end
                if (obj.m_vectorize==0)
                    order = reshape([1:4;5:8],1,8); % alternate u and v
                    KE(:,:,matId) = KE(order,order,matId);
                end
            end
        end
        %% INTEGRATE BODY FORCE AT QUADRILATERAL ELEMENT
        function fE = integrateBodyForceOverElemQ4(obj)
            [xi_GQ,eta_GQ,wt_GQ]= obj.GaussQuad();
            NCell = cell(1,length(xi_GQ));
            for i = 1:length(xi_GQ)
                [NCell{i},~] = obj.QuadShapeFunction(xi_GQ(i),eta_GQ(i));
            end
            fE = zeros(obj.m_numScenarios,8);
            xNodes = obj.m_hx*[0,1,1,0];
            yNodes = obj.m_hy*[0,0,1,1];

            for scenarioId = 1:obj.m_numScenarios
                rho = obj.m_materials(1).rho;
                bx = rho * obj.m_acceleration(scenarioId, 1); % Body force in x-direction
                by = rho * obj.m_acceleration(scenarioId, 2); % Body force in y-direction

                for g = 1:length(wt_GQ)
                    N = NCell{g};
                    J = obj.Jacobian(xNodes,yNodes,xi_GQ(g),eta_GQ(g));
                    dJ = det(J);
                    fE(scenarioId,:) = fE(scenarioId,:) + wt_GQ(g)*dJ*[N*bx, N*by];
                end
                order = reshape([1:4;5:8],1,8); % alternate u and v
                fE(scenarioId,:) = fE(scenarioId,order);
            end
        end
        %% COMPUTE DEFORMATION
        function obj = computeDeformation(obj)
            obj.m_def = zeros(obj.m_ny+1,obj.m_nx+1,obj.m_numScenarios);
            for scenarioId = 1:obj.m_numScenarios
                for i = 1:obj.m_nx+1
                    for j = 1:obj.m_ny+1
                        if (~obj.m_existingNodes(j,i)),continue;end
                        nodeId = ((i-1)*(obj.m_ny+1) + j);
                        obj.m_def(j,i,scenarioId) = sqrt(obj.m_sol(2*nodeId-1,scenarioId)^2+obj.m_sol(2*nodeId,scenarioId)^2);
                    end
                end
            end
        end
        %% COMPUTE VON-MISSES STRESS AT ELEMENTS
        function obj = computeElemStresses(obj)
            if (obj.m_vectorize==1)
                obj = obj.computeElemStresses_vectorize();
            else
                obj = obj.computeElemStresses_loop();
            end
        end

        function obj = computeElemStresses_loop(obj)
            % Compute stresses at the center of element
            nelx = obj.m_nx;
            nely = obj.m_ny;
            nElements = nelx * nely;
            obj.m_strainTensor = zeros(nElements,2,2,obj.m_numScenarios);
            obj.m_stressTensor = zeros(nElements,2,2,obj.m_numScenarios);
            obj.m_vonMisesElems = zeros( obj.m_ny, obj.m_nx,obj.m_numScenarios);
            xNodes = [0,1,1,0]*obj.m_hx;
            yNodes = [0,0,1,1]*obj.m_hy;
            xi = 0; eta = 0;
            J = obj.Jacobian(xNodes,yNodes,xi,eta);
            [~,gradN] = obj.QuadShapeFunction(xi,eta);
            B = J'\gradN;
            for scenarioId = 1:obj.m_numScenarios
                for elx = 1:nelx
                    for ely = 1:nely
                        if (~obj.m_existingElems(ely,elx)),continue;end
                        elem = ((elx-1)*(nely) + ely);
                        nodes = [((elx-1)*(nely+1) + ely),...
                            ((elx)*(nely+1) + ely),...
                            ((elx)*(nely+1) + ely+1),...
                            ((elx-1)*(nely+1) + ely+1)];

                        uvalue = obj.m_sol(2*nodes-1,scenarioId);
                        vvalue = obj.m_sol(2*nodes,scenarioId);

                        gradu = B*uvalue(:);
                        gradv = B*vvalue(:);
                        ux = gradu(1);
                        uy = gradu(2);
                        vx = gradv(1);
                        vy = gradv(2);
                        matId = obj.m_materialIndices(ely,elx);
                        D = obj.m_D{matId};
                        sxx = D(1,1)*ux + D(1,2)*vy;
                        syy = D(2,1)*ux + D(2,2)*vy;
                        sxy = 2*D(3,3)*(uy+vx)/2;
                        obj.m_strainTensor(elem,:,:,scenarioId) = [ux (uy+vx)/2; (uy+vx)/2 vy];
                        obj.m_stressTensor(elem,:,:,scenarioId) = [sxx sxy; sxy syy];
                        obj.m_vonMisesElems(ely,elx,scenarioId) = sqrt(sxx*sxx + syy*syy ...
                            - sxx*syy + 3*sxy*sxy);
                    end
                end
            end
        end

        function obj = computeElemStresses_vectorize(obj)
            % Compute stresses at the center of elements
            nelx = obj.m_nx;
            nely = obj.m_ny;
            nElements = nelx * nely;

            % Pre-allocate tensors for strain, stress, and von Mises stresses
            obj.m_strainTensor = zeros(nElements, 2, 2, obj.m_numScenarios);
            obj.m_stressTensor = zeros(nElements, 2, 2, obj.m_numScenarios);
            obj.m_vonMisesElems = zeros(obj.m_ny, obj.m_nx, obj.m_numScenarios);

            % Define nodal positions for elements
            xNodes = [0, 1, 1, 0] * obj.m_hx;
            yNodes = [0, 0, 1, 1] * obj.m_hy;
            xi = 0; eta = 0;

            % Compute Jacobian and shape function gradient
            J = obj.Jacobian(xNodes, yNodes, xi, eta);
            [~, gradN] = obj.QuadShapeFunction(xi, eta);
            B = J' \ gradN;  % Strain-displacement matrix

            % Get valid element locations
            [elemRows, elemCols] = find(obj.m_existingElems);
            elements = (elemCols - 1) * nely + elemRows;  % Element IDs

            % Vectorize over scenarios
            for scenarioId = 1:obj.m_numScenarios
                % Collect node IDs for each element
                nodes = [(elemCols - 1) * (nely + 1) + elemRows, ...
                    (elemCols) * (nely + 1) + elemRows, ...
                    (elemCols) * (nely + 1) + elemRows + 1, ...
                    (elemCols - 1) * (nely + 1) + elemRows + 1];

                % Iterate over elements since matrix multiplications involve small dimensions
                for idx = 1:numel(elements)
                    elx = elemCols(idx);
                    ely = elemRows(idx);
                    elem = elements(idx);

                    % Get displacement values for the element (u and v components)
                    uvalue = obj.m_sol(2 * nodes(idx, :) - 1, scenarioId);
                    vvalue = obj.m_sol(2 * nodes(idx, :) , scenarioId);

                    % Calculate gradients
                    gradu = B * uvalue(:);
                    gradv = B * vvalue(:);

                    ux = gradu(1);
                    uy = gradu(2);
                    vx = gradv(1);
                    vy = gradv(2);

                    % Material properties
                    matId = obj.m_materialIndices(ely, elx);
                    D = obj.m_D{matId};

                    % Stress components
                    sxx = D(1, 1) * ux + D(1, 2) * vy;
                    syy = D(2, 1) * ux + D(2, 2) * vy;
                    sxy = 2 * D(3, 3) * (uy + vx) / 2;

                    % Store strain and stress tensors
                    obj.m_strainTensor(elem, :, :, scenarioId) = [ux (uy + vx) / 2; (uy + vx) / 2 vy];
                    obj.m_stressTensor(elem, :, :, scenarioId) = [sxx sxy; sxy syy];

                    % Compute von Mises stress
                    obj.m_vonMisesElems(ely, elx, scenarioId) = sqrt(sxx^2 + syy^2 - sxx * syy + 3 * sxy^2);
                end
            end
        end

        %% COMPUTE PRINCIPAL STRESS AT ELEMENTS
        function obj = computePrincipalStress(obj)
            nelx = obj.m_nx;
            nely = obj.m_ny;
            obj.m_principalStress.x =  zeros(nely,nelx,obj.m_numScenarios);
            obj.m_principalStress.y = zeros(nely,nelx,obj.m_numScenarios);
            obj.m_principalStress.tension.value =  zeros(nely,nelx,obj.m_numScenarios);
            obj.m_principalStress.tension.theta =  zeros(nely,nelx,obj.m_numScenarios);
            obj.m_principalStress.tension.u = zeros(nely,nelx,obj.m_numScenarios);
            obj.m_principalStress.tension.v =  zeros(nely,nelx,obj.m_numScenarios);
            obj.m_principalStress.compression.value = zeros(nely,nelx,obj.m_numScenarios);
            obj.m_principalStress.compression.theta = zeros(nely,nelx,obj.m_numScenarios);
            obj.m_principalStress.compression.u = zeros(nely,nelx,obj.m_numScenarios);
            obj.m_principalStress.compression.v =  zeros(nely,nelx,obj.m_numScenarios);

            for scenarioId = 1:obj.m_numScenarios
                for elx = 1:nelx
                    for ely = 1:nely
                        if (~obj.m_existingElems(ely,elx)),continue;end
                        obj.m_principalStress.x(ely,elx,scenarioId) = elx;
                        obj.m_principalStress.y(ely,elx,scenarioId) = ely;
                        elem = ((elx-1)*(nely) + ely);
                        sxx = obj.m_stressTensor(elem,1,1,scenarioId);
                        syy = obj.m_stressTensor(elem,2,2,scenarioId);
                        sxy = obj.m_stressTensor(elem,1,2,scenarioId);


                        stress = [sxx sxy; sxy syy];
                        [directions,values] = eig(stress);

                        val_c = values(1,1); % compression
                        direction_c = directions(:,1);

                        val_t = values(2,2); % tension
                        direction_t = directions(:,2);

                        if (abs(val_t) > abs(val_c)) % tension
                            obj.m_principalStress.tension.value(ely,elx,scenarioId) = val_t;
                            obj.m_principalStress.tension.theta(ely,elx) = atan(direction_t(2)/direction_t(1));
                            obj.m_principalStress.tension.u(ely,elx,scenarioId) = direction_t(1);
                            obj.m_principalStress.tension.v(ely,elx,scenarioId) = direction_t(2);
                        else
                            obj.m_principalStress.compression.value(ely,elx,scenarioId) = val_c;
                            obj.m_principalStress.compression.theta(ely,elx) = atan(direction_c(2)/direction_c(1));
                            obj.m_principalStress.compression.u(ely,elx,scenarioId) = direction_c(1);
                            obj.m_principalStress.compression.v(ely,elx,scenarioId) = direction_c(2);
                        end
                    end
                end
            end
        end
        %% COMPLIANCE
        function cx = computeCompliance(obj)
            U = obj.m_sol;
            f = full(obj.m_f); % force vector
            cx = zeros(1,obj.m_numScenarios);
            for scenarioId = 1:obj.m_numScenarios
                cx(scenarioId) = f(:,scenarioId)' * U(:,scenarioId);
            end
        end
    end
end