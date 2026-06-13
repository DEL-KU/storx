classdef triMesher < brep2d
    properties(GetAccess = 'public', SetAccess = 'public')
        m_elementOrder;
        m_mesh;
        m_elemSize;
        m_nodesPerEdge;
        m_nodesPerElement;
        m_numElemsDesired;
        m_numNodes;
        m_numElems;
        m_pseudoDensity; % for topology optimization
        m_boundaryNodes;
    end
    methods
        function obj = triMesher(brep,nElements,elementOrder)
            % brep can be a brep structure or a file name
            obj = obj@brep2d(brep); % call superclass
            if (nargin == 2)
                elementOrder = 'Quadratic';
            end

            if (nElements > 0)
                obj.m_numElemsDesired = nElements;
                obj = obj.initializeMesher(nElements,elementOrder);
            end
        end
        
        function obj = initializeMesher(obj,nElements,elementOrder)
            % Creates a mesh using Matlab's initmesh
            % Then, depending on the shape function order, additional nodes
            % are inserted.
            % obj.m_elemSize = 1.8*abs(sqrt(obj.m_area/(nElements)));
            % % once again convert using a default of 1 arc segment
            % obj = obj.convertBreptoPdeGeom();
            % [p,e,t] = initmesh(obj.m_pdeGeom,'Hmax',obj.m_elemSize,'MesherVersion', 'R2013a');
            % Compute element size
            obj.m_elemSize = 1.8 * sqrt(obj.m_area / nElements);

            % Generate mesh with modern PDE Toolbox
            mesh = generateMesh(obj.m_pdeModel, ...
                'Hmax', obj.m_elemSize, ...
                'GeometricOrder', 'linear');  % or 'quadratic' if you want higher-order

            % Extract (p,e,t) from mesh
            [p,e,t] = meshToPet(mesh);

            %change boundary number back to original brep boundary numbers
            e(5,:) = obj.m_brepSegsMapping(e(5,:));
            % Force the domain numbering to be reverse ordered
            % This way we always know that domain 0 is to the right
            t1 = (e(6,:) < e(7,:)); % These are the elements we need to flip
            tmp = e([1,3,6],t1);
            e([1,3,6],t1) = e([2,4,7],t1);
            e([2,4,7],t1) = tmp;
            obj.m_mesh = struct('p',p,'e',e,'t',t);% linear
            obj.m_elementOrder = elementOrder;
            obj = obj.findDOFs();
        end

        function obj = findDOFs(obj)
            elementOrder = obj.m_elementOrder;
            if (strcmp(elementOrder,'Linear')) % Linear
                obj.m_nodesPerElement = 3;
                obj.m_nodesPerEdge = 2;
            elseif (strcmp(elementOrder,'Quadratic')) % quadratic
                p = obj.m_mesh.p;
                e = obj.m_mesh.e;
                t = obj.m_mesh.t;
                obj.m_nodesPerElement = 6;
                obj.m_nodesPerEdge = 3;
                mesh = struct('p',p,'e',e,'t',t); % temporary
                % add points on all edges
                np = size(mesh.p,2);
                nt = size(mesh.t,2);
                pairs = [mesh.t([1,2],:), mesh.t([2,3],:), mesh.t([3,1],:)];
                toFlip = pairs(1,:) > pairs(2,:);
                pairs(:,toFlip) = flipud(pairs(:,toFlip));
                [pairs,~,map] = unique(pairs','rows');
                map = reshape(map,nt,3)';
                t_xmesh = [mesh.t(1:3,:); map+np];
                newPts = 0.5*(mesh.p(:,pairs(:,1)) + mesh.p(:,pairs(:,2)));
                [dMin,closestPts] = obj.distOfPointsToBrep(newPts);
                h = obj.m_elemSize;
                newPts(1,:) = (dMin < 0.05*h).*closestPts(1,:) + (dMin >= 0.05*h).*newPts(1,:);
                newPts(2,:) = (dMin < 0.05*h).*closestPts(2,:) + (dMin >= 0.05*h).*newPts(2,:);
                p_xmesh = [mesh.p, newPts];
                e_xmesh = mesh.e;
                e = mesh.e([1,2],:);
                toFlip = e(1,:) > e(2,:);
                e(:,toFlip) = flipud(e(:,toFlip));
                [~,map] = ismember(e',pairs,'rows');
                e_xmesh(3,:) = map + np;
                obj.m_mesh = struct('p',p_xmesh,'e',e_xmesh,'t',t_xmesh);
            end
            obj.m_numNodes = size(obj.m_mesh.p,2);
            obj.m_numElems = size(obj.m_mesh.t,2);
            obj.m_pseudoDensity = ones(obj.m_numElems,1);
            bndryNodes = obj.m_mesh.e(1:2,:);
            bndryNodes = bndryNodes(:);
            obj.m_boundaryNodes = unique(bndryNodes);
        end

        function obj = resetBrepAndMesh(obj,brep)
            % Resets a brep (useful during shape optimization) and
            % recreates the mesh.
            obj = obj.resetBrep(brep);
            obj = obj.initializeMesher(obj.m_numElemsDesired,obj.m_elementOrder);
        end

        function obj = overwritePoints(obj,p)
            % overwrites the mesh pts,  no checks are made
            obj.m_mesh.p = p; % assumes mesh topology does not change
        end

        function obj = snapNodesToBRep(obj, tol)
            % Snaps boundary nodes onto the nearest point on the B-rep
            % while leaving connectivity (t, e) unchanged.
            % tol: only snap nodes whose distance to the B-rep is less than
            %      tol (default: inf, i.e. snap all boundary nodes).
            if nargin < 2
                tol = inf;
            end
            bndryNodes = obj.m_boundaryNodes;
            pts = obj.m_mesh.p(:, bndryNodes);
            [dMin, closestPts] = obj.distOfPointsToBrep(pts);
            toSnap = dMin < tol;
            obj.m_mesh.p(1, bndryNodes(toSnap)) = closestPts(1, toSnap);
            obj.m_mesh.p(2, bndryNodes(toSnap)) = closestPts(2, toSnap);
        end

        function obj = findEdges(obj)
            % given the points and triangles, find the edges on the
            % boundary
            mesh = obj.m_mesh;
            brep = obj.m_brep;
            v = brep.vertices;
            pairs = [mesh.t([1,2],:), mesh.t([2,3],:), mesh.t([3,1],:)];
            toFlip = pairs(1,:) > pairs(2,:);% make sure the first node is smaller
            pairs(:,toFlip) = flipud(pairs(:,toFlip));
            [~,IF] = unique(pairs','rows','first'); % find the unique pairs of nodes
            [~,IL] = unique(pairs','rows','last'); % find the unique pairs of nodes
            index = intersect(IF,IL);
            pairs = pairs(:,index);
            % the unique pairs correspond to edges
            % We need to find the boundary segment the edges fall on
            nMeshEdges = size(pairs,2);
            xStart = mesh.p(1,pairs(1,:));yStart = mesh.p(2,pairs(1,:));
            xEnd = mesh.p(1,pairs(2,:));yEnd = mesh.p(2,pairs(2,:));
            % find the midpoint of each mesh edge
            pts(1,:) = (xStart + xEnd)/2;
            pts(2,:) = (yStart + yEnd)/2;
            % find the distance to each brep segment
            dist = 1e12*ones(obj.m_numBndrySegs,nMeshEdges);
            brep = obj.m_brep;
            for seg = 1:obj.m_numBndrySegs
                breptype = brep.segments(1,seg);
                vs = brep.segments(2,seg);
                ve = brep.segments(3,seg);
                if (breptype == 1) % line
                    dist(seg,:) = obj.distOfPointsToLineSegment(pts,v(:,vs),v(:,ve));
                elseif (breptype == 2) % arc
                    dist(seg,:) = obj.distOfPointsToLineSegment(pts,v(:,vs),v(:,ve));
                    %disp('Distance to arc not implemented');
                end
            end
            obj.m_nodesPerEdge = 2;
            [minDist,bndryIndex] = min(dist);
            if (max(minDist) > 1e-2*obj.m_modelScale)
                disp('Warning: mesh and geometry seem to be inconsistent')
            end
            obj.m_mesh.e = zeros(7,nMeshEdges);
            obj.m_mesh.e(1:2,:) = pairs;
            obj.m_mesh.e(5,:) = bndryIndex;
            obj.m_pseudoDensity = ones(obj.m_numElems,1);
        end

        function obj = setMesh(obj,p,t)
            % overwrites the mesh,  no checks are made
            obj.m_mesh.p = p;
            obj.m_mesh.t = t;
            obj.m_numNodes = size(p,2);
            obj.m_numElems = size(t,2);
            obj = obj.findEdges();
        end
        function obj = readMesh(obj,fileName)
            [fid,message] = fopen(fileName,'r'); % open the file for reading
            if ( 0 > fid) % check to see if the file opened correctly...
                disp(message); % ...and display any error messages
                error('ERROR: failed to open file'); % halt the programme
            end
            nNodes = fscanf(fid,' %d ',1);
            xy = fscanf(fid,' %f ',[2, nNodes]);
            obj.m_numNodes = nNodes;
            obj.m_mesh.p = zeros(2,nNodes);
            obj.m_mesh.p(1:2,:) = xy(1:2,:);
            nElems = fscanf(fid,' %d ',1);
            obj.m_numElems = nElems;
            elems = fscanf(fid,' %f ',[3, nElems]);
            obj.m_mesh.t = ones(4,nElems);
            obj.m_mesh.t(1:3,:) = elems(1:3,:);
            fclose(fid);
            obj = obj.findEdges();
            obj = obj.findDOFs();
        end
        function obj = readOffFileMesh(obj,fileName)
            [fid,message] = fopen(fileName,'r'); % open the file for reading
            if ( 0 > fid) % check to see if the file opened correctly...
                disp(message); % ...and display any error messages
                error('ERROR: failed to open file'); % halt the programme
            end
            str = fgets(fid);   % -1 if eof
            if ~strcmp(str(1:3), 'OFF')
                error('The file is not a valid OFF one.');
            end
            temp = fscanf(fid,' %d ',3);
            nNodes = temp(1);
            nElems = temp(2);
            xy = fscanf(fid,' %f ',[3, nNodes]);
            obj.m_numNodes = nNodes;
            obj.m_mesh.p = zeros(2,nNodes);
            obj.m_mesh.p(1:2,:) = xy(1:2,:);
            obj.m_numElems = nElems;
            elems = fscanf(fid,' %f ',[4, nElems]);
            obj.m_mesh.t = ones(4,nElems);
            obj.m_mesh.t(1:3,:) = elems(2:4,:)+1;%% offset
            fclose(fid);
            obj = obj.findEdges();
            obj = obj.findDOFs();
        end
        function obj = setPseudoDensity(obj,elems,value)
            % Pseudo-density is useful for turning off elements; useful
            % during topology optimization
            obj.m_pseudoDensity(:) = 1;
            obj.m_pseudoDensity(elems) = value;
        end
        function printMesh(obj,fileName)
            % For exporting a mesh
            fptr = fopen(fileName,'w');
            fprintf(fptr,'%d \n',obj.m_numNodes);
            fprintf(fptr,'%f %f \n',obj.m_mesh.p);
            fprintf(fptr,'%d \n',obj.m_numElems);
            fprintf(fptr,'%d %d %d \n',obj.m_mesh.t(1:3,:));
            fclose(fptr);
        end
        function Area = computeMeshArea(obj)
            Area = obj.trimeshArea(obj.m_mesh.p,obj.m_mesh.t);
        end
        function area = computeArea(obj)
            area = obj.computeMeshArea();

            %area = brepArea(obj); %% bug in code
        end
        function Area = trimeshArea(~,p,t)
            % Compute the area of a mesh
            nTriangles = size(t,2);
            Area = 0;
            for elem = 1:nTriangles
                nodes = t(1:3,elem)';
                xNodes = p(1,nodes);
                yNodes = p(2,nodes);
                invJ = [(-yNodes(1)+yNodes(3)) (-yNodes(2)+yNodes(1)); ...
                    (-xNodes(3)+xNodes(1)) (-xNodes(1)+xNodes(2))];
                dJ = invJ(1,1)*invJ(2,2)-invJ(1,2)*invJ(2,1);
                if (dJ < 0)
                    disp('Determinant error in TriMesher');
                end
                Area = Area + dJ/2;
            end
        end
        function obj = setMeshPoints(obj,nodes,p)
            obj.m_mesh.p(:,nodes) = p;
        end
        function  [element,xi,eta] = invert(obj,x,y)
            % Given (x,y) find the element is belongs to, and (xi,eta)
            % not efficient code
            t = obj.m_mesh.t;
            p = obj.m_mesh.p;
            nTriangles = obj.m_numElems;
            element = -1; xi = 0; eta = 0;
            for elem = 1:nTriangles
                nodes = t(1:3,elem)';
                xNodes = p(1,nodes);
                yNodes = p(2,nodes);
                J = [xNodes(2)-xNodes(1)  xNodes(3)-xNodes(1); ...
                    yNodes(2)-yNodes(1)  yNodes(3)-yNodes(1)];
                sol = J\[x-xNodes(1);y-yNodes(1)];
                xi = sol(1);
                eta = sol(2);
                if ((xi >= 0) && (eta >= 0) && (xi+eta<= 1))
                    element = elem;
                    return;
                end
            end
        end
        function  [N,gradN] = triShapeFunction(obj,xi,eta)
            % The shape functions supported
            if (strcmp(obj.m_elementOrder,'Linear')) % linear triangle
                N = [1-xi-eta
                    xi
                    eta];
                gradN = [-1 1 0; -1 0 1];
            elseif (strcmp(obj.m_elementOrder,'Quadratic')) % quadratic triangle
                lambda = 1-xi-eta;
                N = [lambda*(2*lambda-1)
                    xi*(2*xi-1)
                    eta*(2*eta-1)
                    4*xi*lambda
                    4*xi*eta
                    4*eta*lambda];
                gradN = [-3+4*xi+4*eta 4*xi-1 0 4-8*xi-4*eta 4*eta -4*eta;
                    -3+4*xi+4*eta 0 4*eta-1 -4*xi 4*xi 4-4*xi-8*eta];
            end
        end
        function  [fx,fy] = findFieldGradientTri(obj,elem,fieldAtNodes,xi,eta)
            % given nodal values of f, find df/dx and df/dy at any xi, eta
            [~,gradN] = obj.triShapeFunction(xi,eta);
            [J] = obj.Jacobian(elem);
            fieldGradient = J'\(gradN*fieldAtNodes(:));
            fx = fieldGradient(1);
            fy = fieldGradient(2);
        end

        function J = Jacobian(obj,elem)
            t = obj.m_mesh.t;
            p = obj.m_mesh.p;
            nodes = t(1:3,elem)';
            xNodes = p(1,nodes);
            yNodes = p(2,nodes);
            J = [xNodes(2)-xNodes(1)  xNodes(3)-xNodes(1); ...
                yNodes(2)-yNodes(1)  yNodes(3)-yNodes(1)];
        end
        function [N,gradN] = edgeShapeFunction(obj,xi)
            % The edge shape functions supported.
            if (strcmp(obj.m_elementOrder,'Linear'))
                N = [(1-xi)/2;
                    (1+xi)/2];
                gradN = [-1/2;
                    1/2];
            elseif strcmp(obj.m_elementOrder,'Quadratic')
                N = [xi*(xi-1)/2;
                    (1-xi)*(1+xi);
                    xi*(xi+1)/2];
                gradN = [(2*xi-1)/2;
                    (-2*xi);
                    (2*xi+1)/2];
            end
        end
        function plotMesh(obj)
            plt = PlotId;
            figure(plt.mesh); % Create a new figure for the mesh
            % Calls Matlab's pdemesh to plot a mesh
            elemsToPlot = obj.m_pseudoDensity > 0;
            mesh = obj.m_mesh;
            t = mesh.t(1:3,elemsToPlot);
            t(4,:) = 0;
            h = pdemesh(mesh.p,mesh.e,t);
            set(h,'Color','k');
            
            obj.adjustFigScale();
            set(gcf, 'Name', 'Mesh'); % Set the figure name
            xlabel('$x$'); % Label x-axis
            ylabel('$y$'); % Label y-axis
            grid on; % Turn on the grid
            axis on;axis tight; view(2);
            hold off;
            drawnow;
            pause(1e-4);
        end
    end
    methods(Static)
        % Static methods do not require access to member variables
        function [xi_GQ, wt_GQ] = GaussQLine(nPoints)
            % Gauss quadrature pts for a line (-1 to 1)
            if (nPoints <= 1)
                xi_GQ = 0.0;
                wt_GQ = 2;
            elseif (nPoints == 2)
                xi_GQ = [-0.577350269189626 0.577350269189626];
                wt_GQ = [1 1];
            elseif (nPoints == 3)
                xi_GQ = [-0.774596669241483 0 0.774596669241483];
                wt_GQ = [ 0.555555555555556 0.88888888888888889 0.5555555555555556];
            elseif (nPoints == 4)
                xi_GQ = [-0.8611363115 -0.3399810435 0.3399810435 0.8611363115];
                wt_GQ = [0.3478548451 0.6521451548 0.6521451548 0.3478548451];
            elseif (nPoints >= 5)
                xi_GQ = [-0.906179846 -0.53846931 0  0.53846931 0.906179846];
                wt_GQ = [0.236926885 0.4786286704 0.5688888888 0.4786286704 0.236926885];
            end
        end
        function [xi_GQ,eta_GQ, wt_GQ] = GaussQTriangle(nPoints)
            if (nPoints == 1)
                xi_GQ = 1/3;
                eta_GQ = 1/3;
                wt_GQ = 1/2;
            elseif (nPoints == 3)
                xi_GQ = [1/2 1/2 0];
                eta_GQ = [1/2 0 1/2];
                wt_GQ = [1/6 1/6 1/6];
            elseif (nPoints == 7)
                xi_GQ = [0.1012865073235
                    0.7974269853531
                    0.1012865073235
                    0.4701420641051
                    0.4701420641051
                    0.0597158717898
                    0.3333333333333];
                eta_GQ = [0.1012865073235
                    0.1012865073235
                    0.7974269853531
                    0.0597158717898
                    0.4701420641051
                    0.4701420641051
                    0.3333333333333];
                wt_GQ = [0.0629695902724
                    0.0629695902724
                    0.0629695902724
                    0.0661970763942
                    0.0661970763942
                    0.0661970763942
                    0.1125000000000];
            end
        end
    end
end
