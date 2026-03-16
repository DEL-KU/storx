%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: This class defines a grid mesher for 2D boundary             %
% representation (brep) geometries. It generates a structured grid mesh     %
% based on the specified number of elements and whether a uniform grid is   %
% desired. The grid mesher computes the necessary parameters such as        %
% element sizes, node coordinates, and identifies existing nodes and        %
% elements within the geometry. The generated mesh can be used for various  %
% applications in shape and topology optimization.                          %
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

classdef gridMesher <  brep2d
    properties(GetAccess = 'public', SetAccess = 'protected')
        % nx, ny : number of elements in x and y
        % n : number of background elements, n = nx * ny
        % Nx, Ny : number of nodes in x and y, Nx = nx + 1
        % N : number of background nodes, N = Nx * Ny
        m_nx; % number of elements in x
        m_ny; % number of elements in y
        m_hx; % size of elements in x
        m_hy; % size of elements in y
        m_ve; % element area
        m_numElems; % number of all elements Nx x Ny
        m_numNodes; % number of all nodes (Nx+1) x (Ny+1)
        m_existingNodes; % binary matrix of size Nx x Ny
        m_nodeCoords; % coordinates of grid nodes of size 2 x N
        m_elemCoords; % coordinates of element center points n x 2
        m_elemSize; % min size of element
        m_existingElems; % existence of elements (ny x nx) image matrix
        m_numExistingElems; % number of existing elements
        m_bndryNodes; % boundary nodes
        m_bndryElems; % boundary elements
        m_q; % quad elements, each column is a quad element
        m_edges; % edges of elements
        m_nodesPerEdge; % number of nodes per edge
        m_isNodeOnBoundary; % flag indicating if a node is on the boundary
        m_uniformGrid = 0; % whether dx=dy for each grid element (needed for Hamilton-Jacobi solver)
    end
    methods
        %% GRID MESHER CONSTRUCTOR
        function obj = gridMesher(brep,numElements,uniformGrid)
            % GRIDMESHER Constructor for creating a grid mesh from a boundary representation (brep)
            %
            % Input Arguments:
            %     brep - boundary representation structure or file name
            %     numElements - desired number of elements in the mesh
            %     uniformGrid - flag indicating if a uniform grid is desired (default is 0)
            %
            % Output Arguments:
            %     obj - instance of the grid mesh object

            if nargin<3,uniformGrid=0;end % Set default for uniformGrid if not provided
            obj = obj@brep2d(brep); % call superclass constructor
            obj.m_uniformGrid = uniformGrid; % Store uniform grid flag
            h_approx = sqrt(obj.m_area/numElements); % Calculate approximate element size assuming perfect squares
            boundingBoxRatio = obj.m_boxSizes(2)/obj.m_boxSizes(1); % Calculate aspect ratio of bounding box
            obj.m_nx = round(obj.m_boxSizes(1)/h_approx); % Determine number of elements in x-direction
            obj.m_ny = round(obj.m_nx * boundingBoxRatio); % Determine number of elements in y-direction
            obj.m_numElems = obj.m_nx * obj.m_ny; % Calculate total number of elements
            obj.m_numNodes = (obj.m_nx+1) * (obj.m_ny+1); % Calculate total number of nodes
            obj.m_nodesPerEdge = 2; % Set number of nodes per edge
        end
        %% GENERATE GRID MESH
        % generate a grid mesh without conforming to the geometry
        % no assumption is made about the geometry
        function obj = generateGrid(obj)
            % Find the bounding box of the brep
            xMin = obj.m_boundingBox(1,1);  xMax = obj.m_boundingBox(1,2);
            yMin = obj.m_boundingBox(2,1);  yMax = obj.m_boundingBox(2,2);
            % Compute element sizes based on the bounding box and grid dimensions
            obj.m_hx = (xMax-xMin)/obj.m_nx-10*eps;
            obj.m_hy = (yMax-yMin)/obj.m_ny-10*eps;
            if (obj.m_uniformGrid==1)
                % If a uniform grid is specified, set hx equal to hy
                obj.m_hy = obj.m_hx;
                obj.m_boxSizes(1) = obj.m_nx*obj.m_hx;
                obj.m_boxSizes(2) = obj.m_ny*obj.m_hy;
                % Update bounding box dimensions
                obj.m_boundingBox(1,2) = obj.m_boundingBox(1,1) + obj.m_boxSizes(1);
                obj.m_boundingBox(2,2) = obj.m_boundingBox(2,1) + obj.m_boxSizes(2);
            end
            % Determine the minimum element size and volume
            obj.m_elemSize = min(obj.m_hx,obj.m_hy);
            obj.m_ve = obj.m_hx*obj.m_hy;
            % Perturb starting point to avoid numerical issues
            xStart = xMin+5*eps;
            yStart = yMin+5*eps;
            % Create a grid of points and find all points inside the geometry
            [X,Y] = meshgrid((0:obj.m_nx)*obj.m_hx +xStart ,(0:obj.m_ny)*obj.m_hy + yStart);
            obj.m_nodeCoords(1,:) = X(:)';   obj.m_nodeCoords(2,:) = Y(:)';
            inNodes = obj.inBrep(obj.m_nodeCoords); % Find points inside the geometry
            inNodes = reshape(inNodes,obj.m_ny+1,obj.m_nx+1);
            % Include nodes outside that are very close to the brep
            outIndex = find(inNodes == 0);
            ptsOut = obj.m_nodeCoords(:,outIndex);
            dMin = obj.distOfPointsToBrep(ptsOut);
            subIndex = dMin < 0.5*obj.m_elemSize; % Identify nodes close to the brep
            goodNodes = outIndex(subIndex);
            inNodes(goodNodes) = 1; % Mark these nodes as inside
            obj.m_existingNodes = zeros(obj.m_ny+1,obj.m_nx+1);
            obj.m_existingNodes(inNodes>0) = 1; % Create a matrix of existing nodes
            %
            % Find existing solid elements, where all 4 adjacent nodes are solid
            obj.m_existingElems = zeros(obj.m_ny,obj.m_nx);
            obj.m_q = zeros(4,obj.m_nx*obj.m_ny);
            existingElemId = 0;
            for ely = 1:obj.m_ny
                for elx = 1:obj.m_nx
                    elemId = (elx-1)*(obj.m_ny) + ely;
                    % Calculate the coordinates of the element's center
                    obj.m_elemCoords(1,elemId) = xStart + (elx-1)*obj.m_hx + 0.5*obj.m_hx;
                    obj.m_elemCoords(2,elemId) = yStart + (ely-1)*obj.m_hy + 0.5*obj.m_hy;
                    % Determine if the element is solid based on its nodes
                    obj.m_existingElems(ely,elx) = min( ...
                        [obj.m_existingNodes(ely,elx), ...
                        obj.m_existingNodes(ely+1,elx), ...
                        obj.m_existingNodes(ely,elx+1), ...
                        obj.m_existingNodes(ely+1,elx+1)]);
                    if (obj.m_existingElems(ely,elx)==0),continue;end
                    existingElemId = existingElemId + 1;
                    % Store the node indices for the existing element
                    nodes = [((elx-1)*(obj.m_ny+1) + ely),...
                        ((elx)*(obj.m_ny+1) + ely),...
                        ((elx)*(obj.m_ny+1) + ely+1),...
                        ((elx-1)*(obj.m_ny+1) + ely+1)];
                    obj.m_q(:,existingElemId) = nodes;
                end
            end

            % Trim the array to the actual number of elements
            if existingElemId > 0
                obj.m_q = obj.m_q(:,1:existingElemId);
            else
                obj.m_q = zeros(4,0);
            end
            obj.m_numExistingElems = sum(obj.m_existingElems,'all'); % Count total existing elements
            obj = findBoundaryElems(obj); % Identify boundary elements
            obj = findBoundaryNodes(obj); % Identify boundary nodes
            obj = obj.findEdges(); % Find edges of the elements
        end
        %% EDGES AND BOUNDARIES
        function obj = findEdges(obj)
            % FINDINDGES Method to identify edges and boundary nodes in the mesh
            %
            % Input:
            %     obj - instance of the class containing mesh data
            %
            % Output:
            %     obj - updated instance with edges and boundary information

            % Create pairs of nodes for edges
            pairs = [obj.m_q([1,2],:), obj.m_q([2,3],:), obj.m_q([3,4],:),obj.m_q([4,1],:)];
            % Ensure the first node in each pair is the smaller one
            toFlip = pairs(1,:) > pairs(2,:);
            pairs(:,toFlip) = flipud(pairs(:,toFlip));
            % Find unique pairs of nodes
            [~,IF,~] = unique(pairs','rows','first');
            [~,IL] = unique(pairs','rows','last');
            index = intersect(IF,IL);
            pairs = pairs(:,index);

            % Get node IDs that are truly on the specified boundary
            bndryNodeIds = find(obj.m_bndryNodes);  % already computed correctly

            % Filter pairs that consist entirely of boundary nodes
            isBoundaryPair = all(ismember(pairs', bndryNodeIds), 2);
            pairs = pairs(:, isBoundaryPair');

            % Now mark these as edges
            obj.m_edges = zeros(5, size(pairs,2));
            obj.m_edges(1:2,:) = pairs;
            elemSize = min(obj.m_hx,obj.m_hy); % Minimum element size
            boundaryPairs = pairs;
            nMeshEdges = size(boundaryPairs,2);
            % Calculate the center point and closest segment for each edge
            for i=1:nMeshEdges
                sNode = obj.m_nodeCoords(:,boundaryPairs(1,i));
                eNode = obj.m_nodeCoords(:,boundaryPairs(2,i));
                edgeCenterPt = 0.5*(sNode+eNode); % Compute the center point of the edge
                [d,~,closestSeg] = obj.distOfPointsToBrep(edgeCenterPt); % Find the closest segment to the edge center
                if (d > 0.49*elemSize), continue; end % Skip if the distance is too large
                obj.m_edges(5,i) = closestSeg; % Store the closest segment index
            end
            % Identify unique boundary nodes
            bndryNodes = obj.m_edges(1:2,:);
            bndryNodes = unique(bndryNodes(:)); % Get unique boundary node IDs
            % Mark nodes that are on the boundary
            obj.m_isNodeOnBoundary = zeros(1,obj.m_numNodes);
            obj.m_isNodeOnBoundary(bndryNodes) = 1; % Set boundary nodes indicator
        end
        %% SYMMETRY
        function obj = applySymmetryInX(obj)
            % APPLYSYMMETRYINX Adjusts existing nodes and elements for symmetry in the X direction
            obj.m_existingNodes = round(0.5*(obj.m_existingNodes+flipud(obj.m_existingNodes)));
            obj.m_existingElems = round(0.5*(obj.m_existingElems+flipud(obj.m_existingElems)));
        end
        function obj = applySymmetryInY(obj)
            % APPLYSYMMETRYINY Adjusts existing nodes and elements for symmetry in the Y direction
            obj.m_existingNodes = round(0.5*(obj.m_existingNodes+fliplr(obj.m_existingNodes)));
            obj.m_existingElems = round(0.5*(obj.m_existingElems+fliplr(obj.m_existingElems)));
        end
        %% FIND BOUNDARY ELEMENTS
        function obj = findBoundaryElems(obj)
            % FINDBOUNDARYELEMS Method to identify boundary elements in a grid
            %
            % Input:
            %     obj - object containing grid information
            %
            % Output:
            %     obj - updated object with boundary elements marked

            % Initialize boundary elements array with the same size as existing elements
            obj.m_bndryElems = zeros(size(obj.m_existingElems));
            % Loop through each element in the grid
            for elx = 1:obj.m_nx
                for ely = 1:obj.m_ny
                    % Skip if the current element does not exist
                    if (~obj.m_existingElems(ely,elx)),continue;end
                    isBndry = false; % Flag to determine if the element is a boundary element
                    % Check if the element is on the edge of the grid
                    if (elx==1 || elx==obj.m_nx || ely==1 || ely==obj.m_ny)
                        isBndry = true;
                    else
                        % Check neighboring elements to determine if current element is a boundary
                        if (elx > 1)
                            if (~obj.m_existingElems(ely,elx-1)),isBndry = true;end
                        end
                        if (ely > 1)
                            if (~obj.m_existingElems(ely-1,elx)),isBndry = true;end
                        end
                        if (elx < obj.m_nx)
                            if (~obj.m_existingElems(ely,elx+1)),isBndry = true;end
                        end
                        if (ely < obj.m_ny)
                            if (~obj.m_existingElems(ely+1,elx)),isBndry = true;end
                        end
                    end
                    % Assign the boundary status to the corresponding element
                    obj.m_bndryElems(ely,elx) = isBndry;
                end
            end
        end
        %% FIND BOUNDARY NODES
        function obj = findBoundaryNodes(obj)
            obj.m_bndryNodes = false(size(obj.m_existingNodes));
            ny = obj.m_ny;  % number of elements in y
            nx = obj.m_nx;  % number of elements in x

            for nodey = 1:ny+1
                for nodex = 1:nx+1
                    if ~obj.m_existingNodes(nodey, nodex)
                        continue;
                    end

                    % On the outer grid
                    if nodey == 1 || nodey == ny+1 || nodex == 1 || nodex == nx+1
                        obj.m_bndryNodes(nodey, nodex) = true;
                        continue;
                    end

                    % Check the four surrounding elements for existence
                    topLeft    = obj.m_existingElems(nodey-1, nodex-1);
                    topRight   = obj.m_existingElems(nodey-1, nodex);
                    bottomLeft = obj.m_existingElems(nodey,   nodex-1);
                    bottomRight= obj.m_existingElems(nodey,   nodex);

                    if ~(topLeft && topRight && bottomLeft && bottomRight)
                        obj.m_bndryNodes(nodey, nodex) = true;
                    end
                end
            end
        end

        function [nodes,points] = findNodesOnEdge(obj,seg)
            % FINDNODESONEDGE Method to find nodes and points on the edge of a segment
            %
            % Input Arguments:
            %     obj - the object containing the bounding box and node information
            %     seg - the segment index to analyze
            %
            % Output Arguments:
            %     nodes - array of node IDs found on the edge
            %     points - array of points found on the edge

            nodes = []; % Initialize nodes array
            points = []; % Initialize points array
            % Calculate the size of each element in the x and y directions
            elemSize = min(obj.m_hx,obj.m_hx); % Determine the minimum element size

            % Loop through each potential node position
            for i = 1:obj.m_nx+1
                for j = 1:obj.m_ny+1
                    % Skip if the node does not exist or is not a boundary node
                    if (~obj.m_existingNodes(j,i)),continue;end
                    if (~obj.m_bndryNodes(j,i)),continue;end

                    % Calculate the point's coordinates based on the bounding box
                    pt = [obj.m_boundingBox(1,1)+(i-1)*obj.m_hx
                        obj.m_boundingBox(2,1)+(j-1)*obj.m_hy];
                    segType = obj.m_brep.segments(1,seg); % Get the segment type

                    % Calculate the distance to the segment based on its type
                    if (segType == 1)
                        [d,~] = obj.distOfPointsToLineSegment(pt,seg);
                    elseif (segType == 2)
                        [d,~] = obj.distOfPointsToArcSegment(pt,seg);
                    else
                        d = elemSize; % Default distance if segment type is unknown
                    end

                    % If the distance is less than half the element size, store the point and node ID
                    if (d < 0.5*elemSize)
                        points = [points pt]; %#ok
                        nodeId = ((i-1)*(obj.m_ny+1) + j);
                        nodes = [nodes;nodeId]; %#ok
                    end
                end
            end
        end
        %% PLOT GRID MESH AND PSEUDO-DENSITY
        function plotMesh(obj)
            % PLOTMESH Method to plot the grid mesh
            % Sets outside values to NaN, so they are ignored in the plot

            plt = PlotId; % Get the plot identifier
            % Reshape node coordinates for mesh grid
            X = reshape(obj.m_nodeCoords(1,:),[obj.m_ny+1,obj.m_nx+1]);
            Y = reshape(obj.m_nodeCoords(2,:),[obj.m_ny+1,obj.m_nx+1]);
            F = obj.m_existingNodes;
            F(F==0) = NaN; % Set non-existing nodes to NaN for plotting
            figure(plt.mesh); % Create a new figure for the mesh
            surf(X,Y,F); % Create a surface plot
            colormap(gray); % Set the colormap to gray
            view(2); % Set the view to 2D
            pbaspect(obj.m_boxSizes); % Set the aspect ratio of the plot
            set(gcf, 'Name', 'Mesh'); % Set the figure name
            xlabel('$x$'); % Label x-axis
            ylabel('$y$'); % Label y-axis
            axis on; % Turn on the axis
            axis([obj.m_boundingBox(1,:) obj.m_boundingBox(2,:)]); % Set axis limits
            grid on; % Turn on the grid
            plotBndryNodes = false; % Initialize existing element ID
            if plotBndryNodes
                hold on; %#ok % Hold the current plot
                % Get boundary node coordinates
                xBndry = obj.m_nodeCoords(1,obj.m_isNodeOnBoundary(:)==1);
                yBndry = obj.m_nodeCoords(2,obj.m_isNodeOnBoundary(:)==1);
                plot(xBndry,yBndry,'r*') % Plot boundary nodes in red
            end
        end

        function plotWireMesh(obj,figId)

            X = reshape(obj.m_nodeCoords(1,:), [obj.m_ny+1, obj.m_nx+1]);
            Y = reshape(obj.m_nodeCoords(2,:), [obj.m_ny+1, obj.m_nx+1]);

            % Optionally mask lines outside the domain using node mask
            N = obj.m_existingNodes;
            N(N==0) = NaN;               % NaN breaks lines in plot
            N = reshape(N, [obj.m_ny+1, obj.m_nx+1]);

            figure(figId);

            hold on
            % horizontal lines
            for i = 1:size(X,1)
                plot(X(i,:).*N(i,:), Y(i,:).*N(i,:), 'b-');
            end
            % vertical lines
            for j = 1:size(X,2)
                plot(X(:,j).*N(:,j), Y(:,j).*N(:,j), 'b-');
            end
            hold off
        end

    end
end