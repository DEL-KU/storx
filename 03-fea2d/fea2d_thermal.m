%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

classdef fea2d_thermal < fea2d
    properties(GetAccess = 'public', SetAccess = 'private')
        m_fluxSegs;     % Segments for flux calculations
        m_internalHeat; % Internal heat source
        m_T; % Temperature Field
        m_k; % Heat conductivity
    end
    methods (Access = public)
        function obj = fea2d_thermal(brep,numElements,materials,...
                vectorize,numScenarios,interpolation,penaltyStruct,uniformGrid)

            % set default values
            if nargin < 4; vectorize = false;end
            if nargin < 5, numScenarios = 1;end
            if nargin < 6, interpolation = 'none';end
            if nargin < 7; penaltyStruct = struct('min',1,'max',1,'inc',0);end
            if nargin < 8; uniformGrid = 0; end
            % construct
            numDOFperNode = 1; % T
            obj = obj@fea2d(brep,numElements,numDOFperNode,materials,...
                interpolation,numScenarios,penaltyStruct,uniformGrid); % call superclass
            obj.m_internalHeat = zeros(numScenarios,1);

            obj.m_vectorize = vectorize;
             if obj.m_vectorize 
                %% PREPARE FINITE ELEMENT ANALYSIS
                nelx = obj.m_nx;
                nely = obj.m_ny;
                nodenrs = reshape(1:(nelx+1)*(nely+1),nely+1,nelx+1);
                % Element DOFs vector (only one DOF per node in thermal analysis)
                obj.m_edofVec = reshape(nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
                % Element DOFs matrix (4 nodes per element, 1 DOF per node)
                obj.m_edofMat = repmat(obj.m_edofVec,1,4) + ...
                    repmat([-1 0 nely+1 nely], nelx*nely, 1);
                % Kronecker products for assembling global stiffness matrix
                obj.m_iK = reshape(kron(obj.m_edofMat,ones(4,1))',16*nelx*nely,1);
                obj.m_jK = reshape(kron(obj.m_edofMat,ones(1,4))',16*nelx*nely,1);
            end
        end

        function obj = computeMaterialPropertiesMatrices(obj)
            % store heat conductivity
            for matId = 1:obj.m_numMaterials
                % isotropic
                obj.m_k{matId} = obj.m_materials(matId).k;
            end
            % compute template stiffness matrix
            obj.m_KE = obj.integrateKOverElemQ4();
        end
        %% ASSEMBLE STIFFNESS MATRIX
        function obj = assembleK(obj)
            % Assemble global stiffness matrix
            % The global stiffness matrix is assembled in a sparse format
            % The global stiffness matrix is symmetric, so we only need to
            % store the upper triangular part of the matrix.
            % The global stiffness matrix is assembled in a vectorized format
            % if obj.m_vectorize , otherwise it is assembled in a
            % non-vectorized format.
            % The global stiffness matrix is assembled using the element
            % stiffness matrix and the element DOFs.
            % The element stiffness matrix is computed using the material
            % properties and the element shape functions.
            nelx = obj.m_nx;
            nely = obj.m_ny;
            if ~obj.m_vectorize 
                nDOF = (nelx+1)*(nely+1);
                nElements = nelx*nely;
                nzmax = 16*nElements;
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
                        dof = nodes;
                        temp = dof(ones(1,4),:);
                        colIndex = reshape(temp',1,4^2);
                        rowIndex = reshape(temp,1,4^2);
                        entries = reshape(KElem',1,4^2);
                        RowTriplets(index:index+4^2-1,1) = rowIndex';
                        ColTriplets(index:index+4^2-1,1) = colIndex';
                        EntryTriplets(index:index+4^2-1,1) = entries';
                        index = index+4^2;
                    end
                end
                obj.m_K = sparse(RowTriplets,ColTriplets,EntryTriplets,nDOF,nDOF);
            else
                nelx = obj.m_nx;        nely = obj.m_ny;
                xPhys = obj.getInterpolationCoefficient(obj.m_design);
                sK = reshape(obj.m_KE(:)*xPhys(:)',16*nelx*nely,1);
                obj.m_K = sparse(obj.m_iK,obj.m_jK,sK);
                obj.m_K  = (obj.m_K +obj.m_K')/2;
            end
        end
        %% BOUNDARY CONDITIONS
        function obj = applyInternalHeat(obj,value,scenarioId)
            if (nargin < 4),scenarioId=1;end
            obj.m_internalHeat(scenarioId) = value;
        end
        function obj = applyFlux(obj,boundaryEdges,flux,scenarioId)
            if (nargin < 4),scenarioId=1;end
            obj.m_BCtype(boundaryEdges,scenarioId) = 2;
            obj.m_BCvalue(boundaryEdges,scenarioId) = flux./obj.m_segLengths(boundaryEdges);
        end

        function obj = assembleBC(obj)
            [xi_GQ, wt_GQ] = obj.GaussQLine();
            N1D = cell(1,length(xi_GQ));
            for i = 1:length(xi_GQ)
                N1D{i} = obj.edgeShapeFunction(xi_GQ(i));
            end

            nDOF = obj.m_numDOFs;
            % Assemble surface force (Neumann data)
            obj.m_fluxSegs = [];
            obj.m_fluxSegs{obj.m_numScenarios} = [];

            obj.m_forcedNodes = [];
            obj.m_forcedNodes{obj.m_numScenarios} = [];
            for scenarioId = 1:obj.m_numScenarios
                fBoundary = zeros(nDOF,1);
                isDirichlet = zeros(nDOF,1);
                dirValue = zeros(nDOF,1);
                isDirichlet(obj.m_fixedDOFs) = 1;
                for geomEdge = 1:size(obj.m_brep.segments,2)
                    typeT = obj.m_BCtype(geomEdge,1,scenarioId);
                    valueT = obj.m_BCvalue(geomEdge,1,scenarioId);
                    if (typeT == 1)
                        boundarySegments = find(obj.m_edges(5,:) == geomEdge);
                        for seg = boundarySegments
                            nodes = obj.m_edges(1:obj.m_nodesPerEdge,seg);
                            nodes = unique(nodes(:));
                            obj.m_fixedNodes = [obj.m_fixedNodes;  nodes];
                            dof = nodes;
                            isDirichlet(dof) = 1;
                            dirValue(dof) =  valueT;
                        end
                    end
                    if (typeT == 2 ) && (abs(valueT) > 0)% flux
                        boundarySegments = find(obj.m_edges(5,:) == geomEdge);
                        for seg = boundarySegments
                            nodes = obj.m_edges(1:obj.m_nodesPerEdge,seg);
                            obj.m_fluxSegs{scenarioId} = [obj.m_fluxSegs{scenarioId} ;geomEdge*ones(numel(nodes),1)];
                            obj.m_forcedNodes{scenarioId} = unique([obj.m_forcedNodes{scenarioId} ;nodes]);
                            dof = nodes;
                            fBoundaryElem = obj.integrateOverBoundary(geomEdge,seg,wt_GQ,N1D,scenarioId);

                            fBoundary(dof) = fBoundary(dof) + fBoundaryElem;
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
        function obj = assembleInternalLoad(obj)
            obj.m_fE = obj.integrateBodyForceOverElemQ4();
            obj.m_fBody = zeros(obj.m_numDOFs,obj.m_numScenarios);
            nelx = obj.m_nx;
            nely = obj.m_ny;
            % Assemble internal loads
            for scenarioId = 1:obj.m_numScenarios
                for elx = 1:nelx
                    for ely = 1:nely
                        nodes = [((elx-1)*(nely+1) + ely),...
                            ((elx)*(nely+1) + ely),...
                            ((elx)*(nely+1) + ely+1),...
                            ((elx-1)*(nely+1) + ely+1)];
                        dof = nodes;
                        obj.m_fBody(dof,scenarioId) = obj.m_fBody(dof,scenarioId) + obj.m_fE(scenarioId,:)';
                    end
                end
            end
        end

        function fBoundaryElem = integrateOverBoundary(obj,geomEdge,seg,wt_GQ,N1D,scenarioId)
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
                % x = xNodes*N; % radius for axisymmetric problems
                f =  obj.m_BCvalue(geomEdge,1,scenarioId); % dof is either 1 (u) or 2 (v)
                % if (strcmp(obj.myClass,'PlaneStrain')) || (strcmp(obj.myClass,'PlaneStress'))
                fBoundaryElem = fBoundaryElem + wt_GQ(g)*(L/2)*N*f;
                % elseif (strcmp(obj.myClass,'AxiSymmetric'))
                % fBoundaryElem = fBoundaryElem + x*wt_GQ(g)*(L/2)*N*f;
                % else
                %     disp('Error in integrateOverBoundary');
                % end
            end
        end
        %% POST-PROCESSING
        function obj = postProcess(obj)
            obj = obj.computeTemperature();
        end
        %% COMPUTE 2D NODAL TEMPERATURE
        function obj = computeTemperature(obj)
            obj.m_T = zeros(obj.m_ny+1,obj.m_nx+1,obj.m_numScenarios);
            for scenarioId = 1:obj.m_numScenarios
                for i = 1:obj.m_nx+1
                    for j = 1:obj.m_ny+1
                        if (~obj.m_existingNodes(j,i)),continue;end
                        nodeId = ((i-1)*(obj.m_ny+1) + j);
                        obj.m_T(j,i,scenarioId) = obj.m_sol(nodeId,scenarioId);
                    end
                end
            end
        end
        %% GET 2D NODAL TEMPERATURE
        function T = getTemperature(obj,scenarioId)
            if nargin == 1
                T = obj.m_T;
            else
                T = obj.m_T(:,:,scenarioId);
            end
        end
        %% PLOTTING
        function [obj,legend_fields,legend_labels] = plotBoundaryCondition(obj)
            plt = PlotId;
            plbc = PlotBC;
            for scenarioId = 1:obj.m_numScenarios
                legend_fields = {};
                legend_labels = {};
                % plot grid mesh, sets outside value to NaN, so they are ignored in the plot
                obj.plotGeometry(plt.loading+scenarioId,0);
                set(gcf, 'Name', strjoin({'Boundary Condition',num2str(scenarioId)},' ') );
                hold on;
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
                X_Fixed = X(obj.m_fixedDOFs);
                Y_Fixed = Y(obj.m_fixedDOFs);
                if ~isempty(obj.m_fixedDOFs)
                    fixed_T = plot(X_Fixed',Y_Fixed',plbc.fixed_T.marker,'MarkerEdgeColor', plbc.fixed_T.color,'MarkerSize',10);
                    hold on;
                    % legend
                    legend_fields = [legend_fields;fixed_T ]; %#ok
                    legend_labels = [legend_labels,'fixed $T$']; %#ok
                end

                % mark all flux nodes
                nodes = obj.m_forcedNodes{scenarioId};
                if ~isempty(nodes)
                    nodes = nodes(1:2:end);
                    X_forced = X(nodes);
                    Y_forced = Y(nodes);
                    Fx = obj.m_f(nodes,scenarioId);
                    Fy = obj.m_f(nodes,scenarioId);
                    scale = 0.1*obj.m_modelScale;
                    normF = sqrt(Fx.^2 + Fy.^2);
                    for i = 1:numel(nodes)
                        seg = obj.m_fluxSegs{scenarioId}(i);
                        normal = obj.normalOfSegment(seg);
                        Fx(i) = normal(1)*Fx(i);
                        Fy(i) = normal(2)*Fy(i);
                    end
                    if (normF > 0)
                        Fx = Fx./normF;
                        Fy = Fy./normF;
                        Fx = Fx(:)';
                        Fy = Fy(:)';
                        % plot(X_forced,Y_forced,'Marker','diamond','MarkerFaceColor',plbc.flux.color);hold on;
                        start = [X_forced; Y_forced];
                        stop = start + scale*[Fx; Fy];
                        obj.drawArrow(start',stop',plbc.flux.color,2);hold on;
                        % legend
                        flux = plot(NaN, NaN,plbc.flux.marker,'MarkerEdgeColor',plbc.flux.color,'MarkerFaceColor', plbc.flux.color);
                        legend_fields = [legend_fields;flux ]; %#ok
                        legend_labels = [legend_labels,'heat flux']; %#ok
                    end
                end
                %%
                % internal heat
                ih = obj.m_internalHeat(scenarioId);
                if abs(ih) > 0
                    % Get axis limits
                    x_limits = xlim;
                    y_limits = ylim;
                    start = obj.findPointinEmptyRegion(x_limits,y_limits,[X_Fixed;Y_Fixed]);
                    scale = 0.1*obj.m_modelScale;
                    nPoints = 8;
                    % Generate angles from 0 to 2*pi
                    theta = linspace(0, 2*pi, nPoints+1);
                    theta(end) = [];
                    % Parametric equations for the circle
                    stop = start + scale*[cos(theta);sin(theta)];
                    start = start + 0.2*scale*[cos(theta);sin(theta)];
                    if ih < 0
                        tmp = stop;
                        stop = start;
                        start = tmp;
                    end
                    obj.drawArrow(start',stop',plbc.internal_heat.color,4);hold on;
                    % legend
                    internal_heat = plot(NaN, NaN,plbc.internal_heat.marker,'MarkerEdgeColor',plbc.internal_heat.color,'MarkerFaceColor', plbc.internal_heat.color);
                    legend_fields = [legend_fields;internal_heat ]; %#ok
                    legend_labels = [legend_labels,'internal heat']; %#ok
                end
                legend(legend_fields,legend_labels, ...
                    'Location', 'northeastoutside');
               
                pbaspect(obj.m_boxSizes);axis on;axis tight;
            end
            
        end
        function obj = plotTemperature(obj,method)
            plt = PlotId;
            cm = ColorMaps;
            if (nargin==1), method = 'SurfInterp'; end
            % plot grid mesh, sets outside value to NaN, so they are ignored in the plot
            X = reshape(obj.m_nodeCoords(1,:),[obj.m_ny+1,obj.m_nx+1]);
            Y = reshape(obj.m_nodeCoords(2,:),[obj.m_ny+1,obj.m_nx+1]);
            for scenarioId = 1:obj.m_numScenarios
                F = obj.m_T(:,:,scenarioId);
                F(obj.m_existingNodes==0) = nan;  F(obj.m_solidNodes==0) = nan;
                figure(plt.temperature+scenarioId); surf(X,Y,F); 
                colormap(cm.temperature); view(2); colorbar;
                pbaspect(obj.m_boxSizes);axis on;
                xlim(obj.m_boundingBox(1,:));
                ylim(obj.m_boundingBox(2,:));
                set(gcf, 'Name', strjoin({'Temperature',num2str(scenarioId)},' ') );

                if (strcmp(method,'VoxelModel') == 1)
                    grid on;
                elseif (strcmp(method,'SurfInterp') == 1)
                    grid off; shading interp;
                else
                    disp(['Method ' method 'for plotting is not implemented!']);
                end
            end
        end
        function obj = printThermalResults(obj)
            for scenarioId = 1:obj.m_numScenarios
                disp([ 'scenario: ' num2str(scenarioId) ...
                    ', min. Temperature: ' num2str(min(obj.m_T(:,:,scenarioId),[],'all')) ...
                    ', max. Temperature: ' num2str(max(obj.m_T(:,:,scenarioId),[],'all'))]);
            end
        end
    end
    methods (Access = private)
        %% INTEGRATE STIFFNESS MATRIX AT QUADRILATERAL ELEMENT
        function KE = integrateKOverElemQ4(obj)
            [xi_GQ,eta_GQ,wt_GQ]= obj.GaussQuad();
            gradNCell = cell(1,length(xi_GQ));
            for i = 1:length(xi_GQ)
                [~,gradNCell{i}] = obj.QuadShapeFunction(xi_GQ(i),eta_GQ(i));
            end
            KE = zeros(4,4);
            xNodes = obj.m_hx*[0,1,1,0];
            yNodes = obj.m_hy*[0,0,1,1];
            for g = 1:length(wt_GQ)
                gradN = gradNCell{g};
                J = obj.Jacobian(xNodes,yNodes,xi_GQ(g),eta_GQ(g));
                dJ = det(J);
                T = J'\gradN;
                B = T;
                KE = KE + wt_GQ(g)*dJ*B'*obj.m_k{1}*B;
            end
            if obj.m_vectorize
                order = [4 1 2 3];
                KE = KE(order,order);
            end
        end
        %% INTEGRATE INTERNAL HEAT AT QUADRILATERAL ELEMENT
        function fE = integrateBodyForceOverElemQ4(obj)
            [xi_GQ,eta_GQ,wt_GQ]= obj.GaussQuad();
            NCell = cell(1,length(xi_GQ));
            for i = 1:length(xi_GQ)
                [NCell{i},~] = obj.QuadShapeFunction(xi_GQ(i),eta_GQ(i));
            end
            fE = zeros(obj.m_numScenarios,4);
            xNodes = obj.m_hx*[0,1,1,0];
            yNodes = obj.m_hy*[0,0,1,1];
            for scenarioId = 1:obj.m_numScenarios
                for g = 1:length(wt_GQ)
                    N = NCell{g};
                    J = obj.Jacobian(xNodes,yNodes,xi_GQ(g),eta_GQ(g));
                    dJ = det(J);
                    fE(scenarioId,:) = fE(scenarioId,:) + wt_GQ(g)*dJ*N*obj.m_internalHeat(scenarioId);
                end
            end
        end
    end
end