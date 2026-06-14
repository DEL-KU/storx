%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                              %
% This class implements the finite element analysis for 2D elasticity       %
% problems based on triangular elements. It inherits from the triMesher     %
% class and provides methods to compute the elasticity matrix, stiffness    %
% matrix, and perform finite element analysis for elasticity problems.      %                                 %
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

classdef triFEA2d_elasticity < triMesher
    properties(GetAccess = 'public', SetAccess = 'public')
        m_numScenarios = 1;
        m_BCtype;
        m_BCvalue;
        m_DOFPerNode;
        m_NumDOF;
        m_DOFPerElem;
        m_Class;
        m_E;
        m_Nu;
        m_D;
        m_BodyForce;
        m_Xi;
        m_Eta;
        m_Wt;
        m_GradN;
        m_K;
        m_F;
        m_C;
        m_fC;
        m_Sol;
        m_U;
        m_V;
        m_def;
        m_StrainElems;
        m_StressElems;
        m_PrincipalStressesElems;
        m_PrincipalStressNodes;
        m_StrainEnergyElems;
        m_VonMisesElems;
        m_VonMisesNodes;
        m_StrainEnergyDensityNodes;
        m_Compliance;
        m_ElementArea;
        m_FixedDOF; % fixed degrees of freedom
        m_FreeDOF; % all free dof
        m_ForcedNodes; % with non-zero forces
        m_FixedNodes;
        m_SolverMethod;  % solver methods
        m_maxDef; % maximum displacement
        m_maxStress; % max von Mises stress in the domain at any instance
        m_N;
        m_principalStress; % principal stress struct
    end
    methods
        % default: use the parent constructor
        function obj = triFEA2d_elasticity(brepFileName,nElements,materials,class,order)
            if (nargin < 2)
                nElements = 500;
                materials.E = 2e9;
                materials.nu = 0.28;
                class = 'PlaneStress';
                order = 'Quadratic';
            elseif (nargin < 3)
                materials.E = 2e9;
                materials.nu = 0.28;
                class = 'PlaneStress';
                order = 'Quadratic';
            elseif (nargin < 4)
                class = 'PlaneStress';
                order = 'Quadratic';
            elseif (nargin < 5)
                order = 'Quadratic';
            end
            obj = obj@triMesher(brepFileName,nElements,order); % call superclass
            % set the default boundary Conditions
            obj.m_BCtype = zeros(obj.m_numBndrySegs,2);% u v
            obj.m_BCvalue = zeros(obj.m_numBndrySegs,2);% u v
            obj.m_DOFPerNode = 2; % (u, v)
            obj.m_DOFPerElem = obj.m_DOFPerNode*obj.m_nodesPerElement;
            obj.m_NumDOF = obj.m_DOFPerNode*obj.m_numNodes;
            obj.m_E = materials.E;
            obj.m_Nu = materials.nu;
            obj.m_BodyForce = [0 0];
            obj.m_Class = class;

            obj.m_SolverMethod = 2;
            obj = obj.computeDMatrix();

            obj.m_principalStress = struct('x',[],'y',[], ...
                'tension',[],'compression',[]);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = preProcess(obj)
            obj = obj.assembleK();
            obj = obj.assembleBC();
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj,success] = solve(obj)
            [obj,success] = obj.solveLinearSystem();
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = postProcess(obj,principalStressesFlag)
            if nargin == 1, principalStressesFlag = false;end
            obj.m_def = sqrt(obj.m_U.^2 + obj.m_V.^2);
            obj.m_maxDef = max(sqrt(obj.m_U.^2 + obj.m_V.^2));
            obj.m_Compliance = (obj.m_Sol)'*obj.m_F;
            obj = obj.computeStresses();
            obj.m_maxStress = max(obj.m_VonMisesElems);

            if principalStressesFlag,obj = obj.computePrincipalStress();end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = readMesh(obj,fileName)
            obj = obj.readMesh@TriMesher(fileName);
            obj.m_NumDOF = obj.m_DOFPerNode*obj.m_numNodes;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = readOffFileMesh(obj,fileName)
            obj = obj.readOffFileMesh@TriMesher(fileName);
            obj.m_NumDOF = obj.m_DOFPerNode*obj.m_numNodes;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = resetBrepAndSolve(obj,brep)
            % Useful during order optimization
            obj = obj.resetBrepAndMesh(brep);
            obj.m_NumDOF = obj.m_DOFPerNode*obj.m_numNodes;
            obj = obj.solveLinearElasticityProblem(); % Solve Primary FEA problem
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = setSolverMethod(obj,method)
            obj.m_SolverMethod = method;
        end
        function obj = setYoungsModulus(obj,E)
            obj.m_E = E;
            obj = obj.computeDMatrix();
        end
        function obj = setPoissonsRatio(obj,nu)
            obj.m_Nu = nu;
            obj = obj.computeDMatrix();
        end

        function obj = computeDMatrix(obj)
            E = obj.m_E;
            nu = obj.m_Nu;
            if (strcmp(obj.m_Class,'PlaneStrain'))
                obj.m_D = E/((1+nu)*(1-2*nu))*[1-nu nu 0; nu 1-nu 0;0 0 (1-2*nu)/2];
            elseif (strcmp(obj.m_Class,'PlaneStress'))
                obj.m_D = E/(1-nu^2)*[1 nu 0; nu 1 0;0 0 (1-nu)/2];
            elseif (strcmp(obj.m_Class,'AxiSymmetric'))
                obj.m_D = E/((1+nu)*(1-2*nu))*[1-nu nu nu 0; nu 1-nu nu 0; nu nu 1-nu 0;0 0 0 (1-2*nu)/2];
            else
                disp('Class undefined ');
                return;
            end
        end

        function obj = applyXForceOnEdge(obj,boundaryEdges,force)
            % convert force into pressure
            obj.m_BCvalue(boundaryEdges,1) = force./obj.m_segLengths(boundaryEdges);
        end
        function obj = applyYForceOnEdge(obj,boundaryEdges,force)
            % convert force into pressure
            obj.m_BCvalue(boundaryEdges,2) = force./obj.m_segLengths(boundaryEdges);
        end
        function obj = fixXOfEdge(obj,boundaryEdges)
            obj.m_BCtype(boundaryEdges,1) = 1;
        end
        function obj = fixYOfEdge(obj,boundaryEdges)
            obj.m_BCtype(boundaryEdges,2) = 1;
        end
        function obj = fixEdge(obj,boundaryEdges)
            obj.m_BCtype(boundaryEdges,1:2) = 1;
        end
        function obj = setForceVector(obj,f)
            obj.m_F = f;
        end
        function obj = applyBodyForce(obj,b)
            obj.m_BodyForce = b;
        end
        function [udof,vdof] = getDof(obj,elem)
            nodes = obj.m_mesh.t(1:obj.m_nodesPerElement,elem)';
            udof = 2*nodes-1;
            vdof = 2*nodes;
        end

        function obj  = assembleK(obj)
            % Assemble K and f (without boundary conditions) for 2D Linear Elasticity
            % mapping is needed when a subset of elements need to be assembled.
            nDOF = obj.m_numNodes*2;
            if (strcmp(obj.m_elementOrder,'Linear'))
                [xi_GQ,eta_GQ,wt_GQ]= obj.GaussQTriangle(1);
            elseif (strcmp(obj.m_elementOrder,'Quadratic'))
                [xi_GQ,eta_GQ, wt_GQ]= obj.GaussQTriangle(3);
            end
            obj.m_Xi = xi_GQ;
            obj.m_Eta = eta_GQ;
            obj.m_Wt = wt_GQ;

            obj.m_DOFPerElem = obj.m_DOFPerNode*obj.m_nodesPerElement;
            nTriangles = obj.m_numElems;
            nzmax = obj.m_DOFPerElem^2*nTriangles;
            RowTriplets = zeros(nzmax,1);
            ColTriplets = zeros(nzmax,1);
            EntryTriplets = zeros(nzmax,1);
            f = zeros(nDOF,1);
            obj.m_N = cell(1,length(xi_GQ));
            obj.m_GradN = cell(1,length(xi_GQ));
            for i = 1:length(xi_GQ)
                [obj.m_N{i},obj.m_GradN{i}] = obj.triShapeFunction(xi_GQ(i),eta_GQ(i));
            end
            index = 1;
            for elem = 1:nTriangles
                nodes = obj.m_mesh.t(1:obj.m_nodesPerElement,elem)';
                [KElem,fElem,AElem] = obj.integrateKOverElem(elem);
                if (obj.m_pseudoDensity(elem) == 0) % during topology optimization
                    KElem = eye(obj.m_DOFPerElem); %avoids singularity in matrix
                    fElem = zeros(obj.m_DOFPerElem,1);
                else
                    KElem = obj.m_pseudoDensity(elem)^3*KElem;
                    fElem = obj.m_pseudoDensity(elem)*fElem;
                end
                obj.m_ElementArea(elem)= AElem;
                dof = [2*nodes-1; 2*nodes];
                dof = reshape(dof,1,obj.m_DOFPerElem);
                temp = dof(ones(1,obj.m_DOFPerElem),:);
                colIndex = reshape(temp',1,obj.m_DOFPerElem^2);
                rowIndex = reshape(temp,1,obj.m_DOFPerElem^2);
                entries = reshape(KElem',1,obj.m_DOFPerElem^2);
                RowTriplets(index:index+obj.m_DOFPerElem^2-1,1) = rowIndex';
                ColTriplets(index:index+obj.m_DOFPerElem^2-1,1) = colIndex';
                EntryTriplets(index:index+obj.m_DOFPerElem^2-1,1) = entries';
                index = index+obj.m_DOFPerElem^2;
                f(dof) = f(dof) + fElem;
            end
            obj.m_K = sparse(RowTriplets,ColTriplets,EntryTriplets,nDOF,nDOF);
            obj.m_F = f;
        end
        function [KElem,fElem] = computeElementStiffness(obj,elem)
            [KElem,fElem] = integrateKOverElem(obj,elem);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [KElem,fElem,AElem] = integrateKOverElem(obj,elem)

            NCell = obj.m_N;
            gradNCell = obj.m_GradN;
            nodes = obj.m_mesh.t(1:obj.m_nodesPerElement,elem)';
            xNodes = obj.m_mesh.p(1,nodes);
            yNodes = obj.m_mesh.p(2,nodes);
            KElem = zeros(obj.m_DOFPerElem,obj.m_DOFPerElem);
            fElem = zeros(obj.m_DOFPerElem,1);
            % Since we are using linear interpolation for the geometry
            % the Jacobian is a constant.
            invJ = [(-yNodes(1)+yNodes(3)) (-yNodes(2)+yNodes(1)); ...
                (-xNodes(3)+xNodes(1)) (-xNodes(1)+xNodes(2))];
            dJ = invJ(1,1)*invJ(2,2)-invJ(1,2)*invJ(2,1);
            AElem = dJ/2;
            invJ = invJ/dJ;
            Z = zeros(1,obj.m_DOFPerElem/2);
            wt_GQ = obj.m_Wt;
            bx = obj.m_BodyForce(1);
            by = obj.m_BodyForce(2);
            for g = 1:length(wt_GQ)
                N = NCell{g};
                gradN = gradNCell{g};
                x = xNodes*N;% note: for axisymmetric, x is the radius

                T1 = invJ(1,:)*gradN;
                T2 = invJ(2,:)*gradN;
                D = obj.m_D;
                if (strcmp(obj.m_Class,'PlaneStrain'))
                    B = [T1 Z; Z T2; T2 T1];
                    KElem = KElem + wt_GQ(g)*dJ*B'*D*B;
                    fElem = fElem + wt_GQ(g)*dJ*[N*bx;N*by];
                elseif (strcmp(obj.m_Class,'PlaneStress'))
                    B = [T1 Z; Z T2; T2 T1];
                    KElem = KElem + wt_GQ(g)*dJ*B'*D*B;
                    fElem = fElem + wt_GQ(g)*dJ*[N*bx;N*by];
                elseif (strcmp(obj.m_Class,'AxiSymmetric'))
                    B = [T1 Z; Z T2;N'/x Z; T2 T1]; % note: for axisymmetric, x is the radius
                    KElem = KElem + x*wt_GQ(g)*dJ*B'*D*B;
                    fElem = fElem + x*wt_GQ(g)*dJ*[N*bx;N*by];
                end
            end
            if (obj.m_DOFPerElem == 6)
                order = reshape([1:3;4:6],1,6);%order = [1 4 2 5 3 6];
            elseif (obj.m_DOFPerElem == 12)
                order = reshape([1:6;7:12],1,12);
            end
            KElem = KElem(order,order);
            fElem = fElem(order);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = assembleBC(obj)
            if (strcmp(obj.m_elementOrder,'Linear'))
                [xi_GQ, wt_GQ] = obj.GaussQLine(1);
                N1D = cell(1,length(xi_GQ));
                for i = 1:length(xi_GQ)
                    N1D{i} = obj.edgeShapeFunction(xi_GQ(i));
                end
            elseif (strcmp(obj.m_elementOrder,'Quadratic'))
                [xi_GQ, wt_GQ] = obj.GaussQLine(3);
                N1D = cell(1,length(xi_GQ));
                for i = 1:length(xi_GQ)
                    temp = obj.edgeShapeFunction(xi_GQ(i));
                    N1D{i} = temp([1 3 2]); % need to change the order
                end
            end
            nDOF = obj.m_NumDOF;
            % Assemble surface force (Neumann data)
            fBoundary = zeros(nDOF,1);
            obj.m_ForcedNodes = [];
            for geomEdge = 1:size(obj.m_brep.segments,2)
                typeu = obj.m_BCtype(geomEdge,1);
                typev = obj.m_BCtype(geomEdge,2);
                valueu = obj.m_BCvalue(geomEdge,1);
                valuev = obj.m_BCvalue(geomEdge,2);
                if (typeu == 0 ) && (abs(valueu) > 0)
                    boundarySegments = find(obj.m_mesh.e(5,:) == geomEdge);
                    for seg = boundarySegments
                        nodes = obj.m_mesh.e(1:obj.m_nodesPerEdge,seg);
                        udof = 2*nodes-1;
                        fBoundaryElem = obj.integrateOverBoundary(geomEdge,seg,wt_GQ,N1D,1);
                        fBoundary(udof) = fBoundary(udof) +fBoundaryElem;
                        obj.m_ForcedNodes = unique([obj.m_ForcedNodes ;nodes]);
                    end
                end
                if (typev == 0 ) && (abs(valuev) > 0)
                    boundarySegments = find(obj.m_mesh.e(5,:) == geomEdge);
                    for seg = boundarySegments
                        nodes = obj.m_mesh.e(1:obj.m_nodesPerEdge,seg);
                        vdof = 2*nodes;
                        fBoundaryElem = obj.integrateOverBoundary(geomEdge,seg,wt_GQ,N1D,2);
                        fBoundary(vdof) = fBoundary(vdof) + fBoundaryElem;
                        obj.m_ForcedNodes = unique([obj.m_ForcedNodes ;nodes]);
                    end
                end
            end
            % Gather Dirichlet boundary conditions
            isDirichlet = zeros(nDOF,1);
            dirValue = zeros(nDOF,1);
            for geomEdge = 1:size(obj.m_brep.segments,2)
                typeu = obj.m_BCtype(geomEdge,1);
                typev = obj.m_BCtype(geomEdge,2);
                valueu = obj.m_BCvalue(geomEdge,1);
                valuev = obj.m_BCvalue(geomEdge,2);
                if (typeu == 1)
                    boundarySegments = find(obj.m_mesh.e(5,:) == geomEdge);
                    for seg = boundarySegments
                        nodes = obj.m_mesh.e(1:obj.m_nodesPerEdge,seg);
                        obj.m_FixedNodes = [obj.m_FixedNodes;nodes];
                        udof = 2*nodes-1;
                        isDirichlet(udof) = 1;
                        dirValue(udof) =  valueu;
                    end
                end
                if (typev == 1)
                    boundarySegments = find(obj.m_mesh.e(5,:) == geomEdge);
                    for seg = boundarySegments
                        nodes = obj.m_mesh.e(1:obj.m_nodesPerEdge,seg);
                        obj.m_FixedNodes = [obj.m_FixedNodes;nodes];
                        vdof = 2*nodes;
                        isDirichlet(vdof) = 1;
                        dirValue(vdof) =  valuev;
                    end
                end
            end
            dirichletDOF = find(isDirichlet(:) == 1);
            obj.m_FixedNodes = unique(obj.m_FixedNodes);
            obj.m_FixedDOF = dirichletDOF;
            nDirichletDOF = length(dirichletDOF);
            C = zeros(nDirichletDOF,nDOF);
            for i = 1:nDirichletDOF
                C(i,dirichletDOF(i)) = 1;
            end
            C = sparse(C);
            fDirichlet = dirValue(dirichletDOF);
            obj.m_F = obj.m_F + fBoundary;
            obj.m_C = C;
            obj.m_fC = fDirichlet;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function fBoundaryElem = integrateOverBoundary(obj,geomEdge,seg,wt_GQ,N1D,dof)
            nodes = obj.m_mesh.e(1:obj.m_nodesPerEdge,seg);
            xNodes = obj.m_mesh.p(1,nodes);
            yNodes = obj.m_mesh.p(2,nodes);
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
                x = xNodes*N; % radius for axisymmetric problems
                f =  obj.m_BCvalue(geomEdge,dof); % dof is either 1 (u) or 2 (v)
                if (strcmp(obj.m_Class,'PlaneStrain')) || (strcmp(obj.m_Class,'PlaneStress'))
                    fBoundaryElem = fBoundaryElem + wt_GQ(g)*(L/2)*N*f;
                elseif (strcmp(obj.m_Class,'AxiSymmetric'))
                    fBoundaryElem = fBoundaryElem + x*wt_GQ(g)*(L/2)*N*f;
                end
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = setConstraintMatrix(obj,C)
            obj.m_C = C;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj] = setDirichletCondition(obj,dirichletConditions)
            obj.m_fC = dirichletConditions(:);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [obj,success] = solveLinearSystem(obj)
            allDOF = 1:obj.m_NumDOF;
            method = obj.m_SolverMethod;
            success = 2;
            if (method == 1) % direct solve with Lagrange multipliers
                scale = max(max(obj.m_K));
                obj.m_C = scale*obj.m_C;
                obj.m_fC = scale*obj.m_fC;
                [nDirichletDOF] = size(obj.m_C,1);
                % Solution of Algebraic Problem
                Z = spalloc(nDirichletDOF,nDirichletDOF,1);
                KBar = sparse([obj.m_K obj.m_C'; obj.m_C Z]);
                fBar = sparse([obj.m_F;  obj.m_fC]);
                soln = KBar \ fBar;
            elseif (method == 2)% direct solve with elimination  of fixed dof
                obj.m_FreeDOF  = setdiff(allDOF,obj.m_FixedDOF);
                % the useful part of the KMatrix
                KTilde = obj.m_K(obj.m_FreeDOF,obj.m_FreeDOF);
                fTilde = obj.m_F(obj.m_FreeDOF);
                % now subtract all the dirichlet values from rhs
                for i = 1:numel(obj.m_FixedDOF)
                    dof = obj.m_FixedDOF(i);
                    if (abs(obj.m_fC(i)) > 0)
                        fTilde = fTilde - obj.m_K(obj.m_FreeDOF,dof)*obj.m_fC(i);
                    end
                end
                soln = zeros(obj.m_NumDOF,1);
                soln(obj.m_FixedDOF) = obj.m_fC; % fixed values
                soln(obj.m_FreeDOF) =  KTilde\fTilde;
            elseif (method == 3)% iterative solve with elimination of fixed dof
                obj.m_FreeDOF  = setdiff(allDOF,obj.m_FixedDOF);
                % the useful part of the KMatrix
                KTilde = obj.m_K(obj.m_FreeDOF,obj.m_FreeDOF);
                fTilde = obj.m_F(obj.m_FreeDOF);
                % now subtract all the dirichlet values from rhs
                for i = 1:numel(obj.m_FixedDOF)
                    dof = obj.m_FixedDOF(i);
                    if (abs(obj.m_fC(i)) > 0)
                        fTilde = fTilde - obj.m_K(obj.m_FreeDOF,dof)*obj.m_fC(i);
                    end
                end
                soln = zeros(obj.m_NumDOF,1);
                soln(obj.m_FixedDOF) = obj.m_fC; % fixed values
                [soln(obj.m_FreeDOF),~, flag] =  obj.CG(KTilde,fTilde);
                if (flag == -1)
                    success = -1;
                    return;
                end
            end
            nDOF = size(obj.m_K,1);
            obj.m_Sol = full(soln(1:nDOF));
            obj.m_U = full(soln(1:2:nDOF));
            obj.m_V = full(soln(2:2:nDOF));
        end
        function J = computeCompliance(obj)
            J = (obj.m_Sol)'*obj.m_F;
        end
        function [sol,nodes,dof] = getElementSolution(obj,elem)
            nodes = obj.m_mesh.t(1:obj.m_nodesPerElement,elem);
            dof = [2*nodes-1; 2*nodes];
            n = obj.m_nodesPerElement;
            order = reshape([1:n;n+1:2*n],1,2*n);
            dof = dof(order);% alternate u,v
            sol = obj.m_Sol(dof);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = computeStresses(obj)
            % Compute stresses at the average of all Gauss points
            nTriangles = obj.m_numElems;
            nNodes = obj.m_numNodes;
            %[~,gradN] = obj.triShapeFunction(1/3,1/3); % gradient at center
            %gradN = obj.m_GradN;
            obj.m_StrainElems = zeros(nTriangles,2,2);
            obj.m_StressElems = zeros(nTriangles,2,2);
            obj.m_StrainEnergyElems = zeros(nTriangles,1);
            obj.m_VonMisesElems = zeros(1,nTriangles);
            obj.m_VonMisesNodes = zeros(1,nNodes);
            obj.m_StrainEnergyDensityNodes = zeros(1,nNodes);
            nElemsConnectedToNode  = zeros(1,nNodes);
            areasOfElemsConnectedToNode = zeros(1,nNodes);
            D = obj.m_D;
            gradNCell = obj.m_GradN;
            for elem = 1:nTriangles
                if (obj.m_pseudoDensity(elem) == 0)
                    continue;
                end
                obj.m_VonMisesElems(elem) = 0;
                nodes = obj.m_mesh.t(1:obj.m_nodesPerElement,elem)';
                nElemsConnectedToNode(nodes) = nElemsConnectedToNode(nodes) + 1;
                areasOfElemsConnectedToNode(nodes) = areasOfElemsConnectedToNode(nodes) + obj.m_ElementArea(elem);
                xNodes = obj.m_mesh.p(1,nodes);
                yNodes = obj.m_mesh.p(2,nodes);
                uvalue = obj.m_U(nodes);
                vvalue = obj.m_V(nodes);
                invJ = [(-yNodes(1)+yNodes(3)) (-yNodes(2)+yNodes(1)); ...
                    (-xNodes(3)+xNodes(1)) (-xNodes(1)+xNodes(2))];
                dJ = invJ(1,1)*invJ(2,2)-invJ(1,2)*invJ(2,1);
                invJ = invJ/dJ;
                for i = 1:numel(gradNCell)
                    B = invJ*gradNCell{i};
                    gradu = B*uvalue;
                    gradv = B*vvalue;
                    ux = gradu(1);
                    uy = gradu(2);
                    vx = gradv(1);
                    vy = gradv(2);
                    strainElems = [ux (uy+vx)/2; (uy+vx)/2 vy];
                    sxx = D(1,1)*strainElems(1,1) + D(1,2)*strainElems(2,2);
                    syy = D(2,1)*strainElems(1,1) + D(2,2)*strainElems(2,2);
                    sxy = 2*D(3,3)*strainElems(1,2);
                    maxVSStress = sqrt(sxx*sxx + syy*syy - sxx*syy + 3*sxy*sxy);
                    if (maxVSStress > obj.m_VonMisesElems(elem))
                        obj.m_StrainElems(elem,:,:) = strainElems;
                        obj.m_StressElems(elem,:,:) = [sxx sxy; sxy syy];
                        obj.m_VonMisesElems(elem) =  maxVSStress;
                    end
                end
                % [ 1 2 4] to avoid double counting shear
                obj.m_StrainEnergyElems(elem,1) = obj.m_StrainElems(elem,[1 2 4])*obj.m_StressElems(elem,[1 2 4])';
                obj.m_VonMisesNodes(nodes) = obj.m_VonMisesNodes(nodes) + obj.m_VonMisesElems(elem);
                obj.m_StrainEnergyDensityNodes(nodes) = obj.m_StrainEnergyDensityNodes(nodes) + obj.m_StrainEnergyElems(elem);
            end
            obj.m_VonMisesNodes =  obj.m_VonMisesNodes./nElemsConnectedToNode;
            obj.m_StrainEnergyDensityNodes =  obj.m_StrainEnergyDensityNodes./areasOfElemsConnectedToNode;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function maxDef = getMaxDeformation(obj)
            maxDef = obj.m_maxDef;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function maxVonMisesStress = getMaxVonMisesStress(obj)
            maxVonMisesStress = max(obj.m_VonMisesElems);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function pNormStress = getPNormStress(obj,pNorm)
            if (nargin == 1)
                pNorm = 4;
            end
            pNormStress = (sum(obj.m_VonMisesElems.^pNorm)).^(1/pNorm);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function obj = computePrincipalStress(obj)
            % Compute principal stresses at each element
            nTriangles = obj.m_numElems;
            nNodes = obj.m_numNodes;
            obj.m_PrincipalStressesElems = zeros(nTriangles,2);
            obj.m_PrincipalStressNodes = zeros(nNodes,2);
            nElemsConnectedToNode = zeros(nNodes,1);
            for elem = 1:nTriangles
                if (obj.m_pseudoDensity(elem) == 0)
                    continue;
                end
                nodes = obj.m_mesh.t(1:obj.m_nodesPerElement,elem)';
                nElemsConnectedToNode(nodes) = nElemsConnectedToNode(nodes) + 1;
                stress = obj.m_StressElems(elem,1:2,1:2);
                stress = reshape(stress,[2 2]);
                obj.m_PrincipalStressesElems(elem,:) = eig(stress)';
                obj.m_PrincipalStressNodes(nodes,1) = obj.m_PrincipalStressNodes(nodes,1) + obj.m_PrincipalStressesElems(elem,1);
                obj.m_PrincipalStressNodes(nodes,2) = obj.m_PrincipalStressNodes(nodes,2) + obj.m_PrincipalStressesElems(elem,2);
            end
            obj.m_PrincipalStressNodes(:,1) =  obj.m_PrincipalStressNodes(:,1)./nElemsConnectedToNode;
            obj.m_PrincipalStressNodes(:,2) =  obj.m_PrincipalStressNodes(:,2)./nElemsConnectedToNode;
        end
        %% COMPUTE PRINCIPAL STRESS AT ELEMENTS
        %% RESULTS
        function obj = printElascticityResults(obj)
            for scenarioId = 1:obj.m_numScenarios
                disp([ 'scenario: ' num2str(scenarioId) ...
                    ', max. Deformation: ' num2str(obj.m_maxDef) ...
                    ', max. vonMises: ' num2str(obj.m_maxStress)]);
            end
        end
        %% PLOTTING

        function plotBoundaryCondition(obj)
            plt = PlotId;
            plbc = PlotBC;
            scenarioId = 1;
            legend_fields = {};
            legend_labels = {};

            obj.plotGeometry(plt.loading+scenarioId,0);
            set(gcf, 'Name', strjoin({'Boundary Condition',num2str(scenarioId)},' ') );hold on;
            % mark all fixed nodes
            index = (rem(obj.m_FixedDOF,2)==1);
            fixedXNodes = (obj.m_FixedDOF(index)-1)/2+1;
            index = (rem(obj.m_FixedDOF,2)==0);
            fixedYNodes = (obj.m_FixedDOF(index))/2;
            % onlyXFixed = setdiff(fixedXNodes,fixedYNodes);
            % onlyYFixed = setdiff(fixedYNodes,fixedXNodes);
            % bothFixed = intersect(fixedXNodes,fixedYNodes);
            fixed_U = plot(obj.m_mesh.p(1,fixedXNodes),obj.m_mesh.p(2,fixedXNodes),plbc.fixed_U.marker,'MarkerEdgeColor', plbc.fixed_U.color);
            % legend
            legend_fields = [legend_fields;fixed_U ]; %#ok
            legend_labels = [legend_labels,'fixed $u$']; %#ok

            fixed_V = plot(obj.m_mesh.p(1,fixedYNodes),obj.m_mesh.p(2,fixedYNodes),plbc.fixed_V.marker,'MarkerEdgeColor', plbc.fixed_V.color);
            % legend
            legend_fields = [legend_fields;fixed_V ]; %#ok
            legend_labels = [legend_labels,'fixed $v$']; %#ok

            % mark all forced nodes
            nodes = obj.m_ForcedNodes;
            Fx = obj.m_F(2*nodes-1);
            Fy = obj.m_F(2*nodes);
            xRange = max(obj.m_mesh.p(1,:))-min(obj.m_mesh.p(1,:));
            yRange = max(obj.m_mesh.p(2,:))-min(obj.m_mesh.p(2,:));
            scale = 0.25*min(xRange,yRange);
            normF = sqrt(Fx.^2 + Fy.^2);
            if (normF > 0)
                Fx = Fx./normF;
                Fy = Fy./normF;
                Fx = Fx(:)';
                Fy = Fy(:)';
                plot(obj.m_mesh.p(1,nodes),obj.m_mesh.p(2,nodes),'Marker','diamond','MarkerFaceColor',plbc.force.color);hold on;
                start = [obj.m_mesh.p(1,nodes); obj.m_mesh.p(2,nodes)];
                stop = start + scale*[Fx; Fy];
                obj.drawArrow(start',stop',plbc.force.color,2);hold on;
                % legend
                force = plot(NaN, NaN,plbc.force.marker,'MarkerEdgeColor',plbc.force.color,'MarkerFaceColor', plbc.force.color); hold on
                legend_fields = [legend_fields;force ]; %#ok
                legend_labels = [legend_labels,'force']; %#ok
            end
            legend(legend_fields,legend_labels, ...
                'Location', 'northeastoutside');

            axis auto;axis equal;
            obj.adjustFigScale();

            hold off;
            drawnow;
            pause(1e-4);
        end
        function printMesh(obj,fileName)
            fptr = fopen(fileName,'w');
            p = zeros(3,obj.m_numNodes);
            p(1:2,:) = obj.m_mesh.p;
            for i = 1:obj.m_numNodes
                if (ismember(i,obj.m_FixedNodes))
                    p(3,i) = 1;
                elseif (ismember(i,obj.m_ForcedNodes))
                    p(3,i) = 2;
                else
                    p(3,i) = 0;
                end
            end
            fprintf(fptr,'%d \n',obj.m_numNodes);
            fprintf(fptr,'%f %f %d\n',p);
            fprintf(fptr,'%d \n',obj.m_numElems);
            if (strcmp(obj.m_elementOrder,'Quadratic'))
                fprintf(fptr,'%d %d %d %d %d %d\n',obj.m_mesh.t(1:6,:));
            else
                fprintf(fptr,'%d %d %d \n',obj.m_mesh.t(1:3,:));
            end
            fclose(fptr);
            hold off;
            drawnow;
            pause(1e-4);
        end
        function plotDeformation(obj)
            plt = PlotId;
            cm = ColorMaps;
            scale = 0.1*obj.m_modelScale/obj.m_maxDef;
            elemsToPlot = obj.m_pseudoDensity > 0;
            mesh = obj.m_mesh;
            scenarioId = 1;
            p = mesh.p + scale*[obj.m_U';obj.m_V'];
            t = mesh.t(1:3,elemsToPlot);
            t(4,:) = 0;
            x = p(1,:);
            y = p(2,:);
            figure(plt.deformation+scenarioId);
            elemsToPlot = obj.m_pseudoDensity > 0;
            trisurf(t(1:3,elemsToPlot)',x,y,obj.m_def,obj.m_def,'facecolor','interp'); hold on
            set(gcf, 'Name', strjoin({'Deformation ',num2str(scenarioId)},' ') );
            pdemesh(mesh.p,mesh.e,mesh.t);
            colormap(cm.deformation); colorbar;
            obj.adjustFigScale();
            axis off;axis tight; grid off; view(2);
            hold off;
            drawnow;
            pause(1e-4);
        end
        function plotVonMisesStress(obj)
            plt = PlotId;
            cm = ColorMaps;
            x = obj.m_mesh.p(1,:);
            y = obj.m_mesh.p(2,:);
            scenarioId = 1;
            vonMises = obj.m_VonMisesNodes;
            elemsToPlot = obj.m_pseudoDensity > 0;
            figure(plt.von_mises+scenarioId);
            trisurf(obj.m_mesh.t(1:3,elemsToPlot)',x,y,vonMises,vonMises,'facecolor','interp');
            colormap(cm.von_mises); hold on
            clim([0 max(vonMises)]);
            set(gcf, 'Name', strjoin({'von Mises Stress ',num2str(scenarioId)},' ') );
            obj.adjustFigScale();
            axis off;axis tight; grid off; view(2);
            pbaspect(obj.m_boxSizes);axis off;
            xlim(obj.m_boundingBox(1,:));
            ylim(obj.m_boundingBox(2,:));
            colorbar;
            hold off;
            drawnow;
            pause(1e-4);
        end
        function plotFirstPrincipalStress(obj)
            x = obj.m_mesh.p(1,:);
            y = obj.m_mesh.p(2,:);
            data = obj.m_PrincipalStressNodes(:,1);
            elemsToPlot = obj.m_pseudoDensity > 0;
            trisurf(obj.m_mesh.t(1:3,elemsToPlot)',x,y,data,data,'facecolor','interp','edgecolor','none');
            axis off;axis tight; grid off; view(2);
            obj.adjustFigScale();
            clim([min(data) max(data)]); colormap('default');
            colorbar('location','eastoutside')
            title('FirstPrincipalStress');
            hold off;
            drawnow;
            pause(1e-4);
        end
        function plotSecondPrincipalStress(obj)
            x = obj.m_mesh.p(1,:);
            y = obj.m_mesh.p(2,:);
            data = obj.m_PrincipalStressNodes(:,2);
            elemsToPlot = obj.m_pseudoDensity > 0;
            trisurf(obj.m_mesh.t(1:3,elemsToPlot)',x,y,data,data,'facecolor','interp','edgecolor','none');
            axis off;axis tight; grid off; view(2);
            obj.adjustFigScale();
            clim([min(data) max(data)]); colormap('default');
            colorbar('location','eastoutside')
            title('SecondPrincipalStress');
            hold off;
            drawnow;
            pause(1e-4);
        end

        function plotStrainEnergyDensity(obj)
            x = obj.m_mesh.p(1,:);
            y = obj.m_mesh.p(2,:);
            energy = obj.m_StrainEnergyDensityNodes;
            elemsToPlot = obj.m_pseudoDensity > 0;
            trisurf(obj.m_mesh.t(1:3,elemsToPlot)',x,y,energy,energy,'facecolor','interp','edgecolor','none');
     
            obj.adjustFigScale();
            axis off;axis tight; grid off; view(2);
            hold off;
            clim([0 max(energy)]); colormap('default');
            colorbar('location','eastoutside')
            title(['StrainEnergyDensity: ' num2str(max(energy))]);
            hold off;
            drawnow;
            pause(1e-4);
        end

        function plotPrincipalStress(obj)
            sigma1 = obj.m_PrincipalStressNodes(:,1);
            sigma2 = obj.m_PrincipalStressNodes(:,2);
            elemsToPlot = obj.m_pseudoDensity > 0;
            tricontour(obj.m_mesh.p',obj.m_mesh.t(1:3,elemsToPlot)',sigma1,20);
            axis off;axis tight; grid off; view(2);
            hold on;
            tricontour(obj.m_mesh.p',obj.m_mesh.t(1:3,elemsToPlot)',sigma2,20);
            obj.adjustFigScale();grid off;
            title('Principal Stresses Contours');
            hold off;
            drawnow;
            pause(1e-4);
        end

    end
end