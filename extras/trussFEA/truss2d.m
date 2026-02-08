classdef truss2d
    properties(GetAccess = 'public', SetAccess = 'private')
        % public read access, but private write access.
        myNodeLocations; % (x,y) node locations,(2, N), N is number of nodes
        myConnectivity; % (startnode, endnode) (2, M), M is number of truss bars
        myNumNodes; % number of Nodes
        myNumTrussBars; % number of truss bars
        myE; % Young's modulus, (1, myNumTrussBars)
        myArea; % cross-sectional area, (1, myNumTrussBars)
        myL; % Length of each bar
        myOrientation; % orientation of each bar, (2,M)
        myDOFFree; % 0 or 1, (1,2N), N is the number of nodes
        myK; % stiffness matrix, (2N, 2N)
        myForceExternal; %force vector at each node, (2N,1)
        mySol; % (2N,1) solution
        myUV; % (2,N) (u,v)  displacements of all nodes
        myDeformations; % deformation of each member
        myDeformationMax; % max deformation among all members
        myStress; % stress in each member
    end
    methods
        function obj = truss2d(nodeXY,connectivity)
            % constructor for a truss system
            % nodeXY: 2 x N, N is the number of nodes
            % connectivity: 2 x M, M is the number of truss bars
            if (size(nodeXY,1) ~=2)
                nodeXY = nodeXY'; % transpose
            end
            if (size(connectivity,1) ~=2)
                connectivity = connectivity'; % transpose
            end
            obj.myNodeLocations = nodeXY;
            obj.myConnectivity = connectivity;
            obj.myNumNodes = size(nodeXY,2);
            obj.myNumTrussBars = size(connectivity,2);
            obj.myE = 2e11*ones(1,obj.myNumTrussBars);
            obj.myArea = 1e-4*ones(1,obj.myNumTrussBars);
            obj.myDOFFree = ones(2*obj.myNumNodes,1);% all free nodes
            obj.myForceExternal = zeros(2*obj.myNumNodes,1);
            obj.myUV = zeros(2,obj.myNumNodes);
            obj.myStress = zeros(1,obj.myNumTrussBars);

            % Find the length and orientation of each bar
            C = obj.myConnectivity;
            startPts = obj.myNodeLocations(:,C(1,:));
            endPts = obj.myNodeLocations(:,C(2,:));
            obj.myL = sqrt((endPts(1,:)-startPts(1,:)).^2 + (endPts(2,:)-startPts(2,:)).^2);
            obj.myOrientation(1,:) = (endPts(1,:)-startPts(1,:))./obj.myL;
            obj.myOrientation(2,:) = (endPts(2,:)-startPts(2,:))./obj.myL;
            obj.myDeformationMax = 0;
        end
        function obj = assignE(obj,E,members)
            % obj = assignE(obj,E,members)
            % assign E to one or more members
            % if members is not give, then assign E to all members
            if (nargin == 2)
                members = 1:obj.myNumTrussBars;
            else
                assert(max(members) <= obj.myNumTrussBars);
                assert(min(members) >=  1);
            end
            obj.myE(members) = E;
        end
        function obj = assignA(obj,A,members)
            %assign Area to one or more members
            % if members = [], then assign A to all members
            if (nargin == 2)
                members = 1:obj.myNumTrussBars;
            else
                assert(max(members) <= obj.myNumTrussBars);
                assert(min(members) >=  1);
            end
            obj.myArea(members) = A;
        end
        function obj = fixXofNodes(obj,nodes)
            %fix the x locations of the nodes
            assert(max(nodes) <= obj.myNumNodes);
            assert(min(nodes) >=  1);
            obj.myDOFFree(2*nodes-1) = 0;
        end
        function obj = fixYofNodes(obj,nodes)
            %fix the y locations of the nodes
            assert(max(nodes) <= obj.myNumNodes);
            assert(min(nodes) >=  1);
            obj.myDOFFree(2*nodes) = 0;
        end
        function obj = applyForce(obj,node,force)
            % apply force at specified nodes
            assert(max(node) <= obj.myNumNodes);
            assert(min(node) >=  1);
            obj.myForceExternal(2*node-1,1) = force(1);% x force
            obj.myForceExternal(2*node,1) = force(2);% y force
        end
        function plotWithNumbering(obj)
            obj.plot(1);
        end

        function obj = assemble(obj)
            obj = obj.assembleElemBased();
        end
        function obj = assembleNodeBased(obj) % Skip
            % Node based assembly of K matrix
            N = obj.myNumNodes; % number of nodes
            C = obj.myConnectivity;
            KMatrix = zeros(2*N,2*N);
            % k = EA/L of each member
            Ke = (obj.myE).*(obj.myArea)./(obj.myL);
            % for each node
            for node = 1:N
                % find all elements connected to it
                % the node can occur in the first or second position
                elems = union(find(C(1,:)== node),find(C(2,:)== node));
                for k = elems % for each elem
                    % find the other node
                    node2 = setdiff(C(:,k),node);
                    c = obj.myOrientation(1,k); % cos(alpha)
                    s = obj.myOrientation(2,k); % sin(alpha)
                    KMatrix(2*node-1,2*node-1) =  KMatrix(2*node-1,2*node-1) +  Ke(k)*c*c;
                    KMatrix(2*node-1,2*node) = KMatrix(2*node-1,2*node) + Ke(k)*c*s;
                    KMatrix(2*node,2*node-1) = KMatrix(2*node,2*node-1) + Ke(k)*c*s;
                    KMatrix(2*node,2*node) = KMatrix(2*node,2*node) + Ke(k)*s*s;

                    KMatrix(2*node-1,2*node2-1) =  -Ke(k)*c*c;
                    KMatrix(2*node-1,2*node2) =  -Ke(k)*c*s;
                    KMatrix(2*node,2*node2-1) =  -Ke(k)*c*s;
                    KMatrix(2*node,2*node2) =  -Ke(k)*s*s;
                end
            end
            obj.myK = KMatrix;
        end
        function obj = assembleElemBased(obj)
            % Simpler and more efficient means of assembly
            N = obj.myNumNodes; % number of nodes
            M = obj.myNumTrussBars;
            C = obj.myConnectivity;
            KMatrix = zeros(2*N,2*N);
            % k = EA/L of each member
            Ke = (obj.myE).*(obj.myArea)./(obj.myL);
            for k = 1:M
                m = C(1,k);
                n = C(2,k);
                c = obj.myOrientation(1,k); % cos(alpha)
                s = obj.myOrientation(2,k); % sin(alpha)
                kElem = Ke(k)*[c^2 c*s ; c*s s^2];
                mDOF = [2*m-1 2*m];
                nDOF = [2*n-1 2*n];
                KMatrix(mDOF,mDOF) = KMatrix(mDOF,mDOF) + kElem;
                KMatrix(nDOF,nDOF) = KMatrix(nDOF,nDOF) + kElem;
                KMatrix(mDOF,nDOF) = KMatrix(mDOF,nDOF) - kElem;
                KMatrix(nDOF,mDOF) = KMatrix(nDOF,mDOF) - kElem;
            end
            obj.myK = KMatrix;
        end
        function obj = solve(obj)
            % solve the truss problem, after the matrix has been assembled
            freeDOF = find(obj.myDOFFree);
            KTilde = obj.myK(freeDOF,freeDOF);
            fTilde = obj.myForceExternal(freeDOF);
            sol = zeros(2*obj.myNumNodes,1);
            sol(freeDOF) = KTilde\fTilde; % solve for only the free dof

            obj.mySol = sol;
            % now compute other useful stuff
            obj.myUV(1,:) = sol(1:2:end);
            obj.myUV(2,:) = sol(2:2:end);
            u = obj.myUV(1,:);
            v = obj.myUV(2,:);
            obj.myDeformations = sqrt(u.^2 + v.^2);
            obj.myDeformationMax = max(obj.myDeformations);
            obj = obj.computeStresses();
        end
        function vol = getVolume(obj)
            vol = sum(obj.myArea.*obj.myL);
        end
        function J = getCompliance(obj)
            J = obj.myForceExternal(:)'*obj.mySol;
        end
        function U = elasticEnergyInBar(obj,k)
            C = obj.myConnectivity;
            u = obj.myUV(1,:);
            v = obj.myUV(2,:);
            m = C(1,k);
            n = C(2,k);
            c = obj.myOrientation(1,k); % cos(alpha)
            s = obj.myOrientation(2,k); % sin(alpha)
            stiffness = obj.myE(k)*obj.myArea(k)/obj.myL(k);
            delta = (u(n)-u(m))*c + (v(n)-v(m))*s;
            U =0.5*stiffness*delta^2;
        end
        function obj = computeStresses(obj)
            % after the truss problem is solved, find the stresses in each element
            M = obj.myNumTrussBars;
            C = obj.myConnectivity;
            Se = (obj.myE)./(obj.myL);
            u = obj.myUV(1,:);
            v = obj.myUV(2,:);
            for k = 1:M % for each member
                m = C(1,k);
                n = C(2,k);
                c = obj.myOrientation(1,k); % cos(alpha)
                s = obj.myOrientation(2,k); % sin(alpha)
                delta = (u(n)-u(m))*c + (v(n)-v(m))*s;
                obj.myStress(k) = Se(k)*delta;
            end
        end

        %% PLOT
        function plot(obj, withNumbering,plotId)
            % Plot the undeformed truss with bar thickness & color mapped to Area
            if nargin < 3, plotId = 1; end
            if nargin < 2,withNumbering = 0;end

            % --- data ---
            A  = obj.myArea(:);                 % Area per bar (Mx1)
            N  = obj.myNumNodes;                % number of nodes
            M  = obj.myNumTrussBars;            % number of bars
            P  = obj.myNodeLocations;           % node XY (2 x N)
            C  = obj.myConnectivity;            % bar connectivity (2 x M)
            L  = obj.myL;                       % bar lengths (1 x M)
            LScale = max(L);

            % --- visuals: linewidth & colormap mapping from Area ---
            lw_min = 1.5;       % min line width
            lw_max = 6.0;       % max line width
            cmap   = cool(256); % use 'turbo' (or 'parula' if older MATLAB)
            Amin = min(A); Amax = max(A);
            if Amax == Amin
                % Avoid division by zero: treat all bars the same
                tA = zeros(M,1);
            else
                tA = (A - Amin) / (Amax - Amin + 1e-8); % normalize to [0,1]
            end

            % --- figure setup ---
            figure(plotId); hold on; axis equal; box on;
            xlabel('$x$','Interpreter','latex'); ylabel('$y$','Interpreter','latex');
            title('Undeformed truss (bar color/width $\propto$ Area)','Interpreter','latex','FontSize',18);

            % --- plot bars: color & thickness by Area ---
            for m = 1:M
                startNode = C(1,m);
                endNode   = C(2,m);
                xTruss = P(1,[startNode endNode]);
                yTruss = P(2,[startNode endNode]);

                % color from colormap
                idx = max(1, min(size(cmap,1), round(1 + tA(m)*(size(cmap,1)-1))));
                col = cmap(idx, :);

                % line width from Area
                lw = lw_min + tA(m)*(lw_max - lw_min);

                plot(xTruss, yTruss, '-', 'Color', col, 'LineWidth', lw);

                if withNumbering
                    xText = 0.4*xTruss(1)+0.6*xTruss(2);
                    yText = 0.4*yTruss(1)+0.6*yTruss(2);
                    text(xText, yText, num2str(m), 'FontAngle','italic', ...
                        'Color',[0.2 0.2 0.2], 'FontSize', 12, 'Interpreter','none');
                end
            end

            % --- nodes (same as your original) ---
            for j = 1:N
                x0 = P(1,j); y0 = P(2,j);
                ufree = obj.myDOFFree(2*j-1);
                vfree = obj.myDOFFree(2*j);
                if (ufree && vfree)               % both dof free
                    h = plot(x0,y0,'ro'); set(h,'LineWidth',2);
                    if withNumbering, text(x0+0.05*LScale,y0,num2str(j),'Color',[1 0 0],'FontSize',12); end
                elseif (ufree && ~vfree)          % one dof free
                    h = plot(x0,y0,'bv'); set(h,'LineWidth',2.0);
                    if withNumbering, text(x0+0.05*LScale,y0,num2str(j),'Color',[0 0 0],'FontSize',12); end
                elseif (~ufree && vfree)          % one dof free
                    h = plot(x0,y0,'b>'); set(h,'LineWidth',2.0);
                    if withNumbering, text(x0+0.05*LScale,y0,num2str(j),'Color',[0 0 0],'FontSize',12); end
                else                               % both dof fixed
                    h = plot(x0,y0,'ko','MarkerFaceColor','k'); set(h,'LineWidth',5.0);
                    if withNumbering, text(x0+0.05*LScale,y0,num2str(j),'Color',[0 0 0],'FontSize',12); end
                end
            end

            % --- external forces (same as your original) ---
            Fx = obj.myForceExternal(1:2:end)';
            Fy = obj.myForceExternal(2:2:end)';
            scale = 0.25;
            hq = quiver(P(1,:), P(2,:), Fx, Fy, scale);
            set(hq,'LineWidth',2.0,'Color','r');

            % --- colorbar keyed to Area ---
            if Amax ~= Amin
            colormap(cmap);
            c = colorbar;
            c.Label.String = 'Bar Area (Normalized)';
            c.Label.Interpreter = 'latex';
            % clim([Amin Amax]); % ensure color limits match Areas

            % optional neat tick labels if ranges are large/small
            % if Amax > 0 && Amax ~= Amin
            %     ticks = linspace(Amin, Amax, 5);
            % 
            %     c.Ticks = ticks;
            %     c.TickLabels = arrayfun(@(x)sprintf('%.3g',x), ticks, 'UniformOutput', false);
            % end
            end
            hold off;
            drawnow
        end

        function plotDeformed(obj,plotId)
            
            % plot the deformed truss
            if (obj.myDeformationMax < eps)
                % solution does not exist
                return;
            end
            figure(plotId);
            obj.plot(0,plotId); hold on
            M = obj.myNumTrussBars;
            P = obj.myNodeLocations; % extract undeformed location
            C = obj.myConnectivity;
            % A template truss element of unit length along x axis
            L = obj.myL;
            LScale = max(L);
            scale = 0.15*LScale/obj.myDeformationMax;
            P = P + scale*obj.myUV;
            for m = 1:M % for all trusses
                startNode = C(1,m);
                endNode = C(2,m);
                xTruss = P(1,[startNode endNode]);
                yTruss = P(2,[startNode endNode]);
                plot(xTruss,yTruss,':k','LineWidth',1.5); hold on;axis('equal');
            end
            v = obj.getVolume();
            title(['$\delta =$ ' num2str(obj.myDeformationMax,'%0.1e') '; $\sigma =$ ' num2str(max(abs(obj.myStress)),'%0.1e') ';  vol= ' num2str(v,'%.2e')],'Interpreter','latex','FontSize',18);
            drawnow
        end


    end
end
