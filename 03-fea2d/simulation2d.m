%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discription:                                                              %       
% This is an abstract class for 2D simulations in the context of shape      %
% and topology optimization. It inherits from the gridMesher class and      %   
% provides methods for pre-processing, solving, and post-processing the     %
% simulation. The class is designed to handle design variables,             %
% materials, and load scenarios in a 2D grid mesh environment.              %
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

classdef (Abstract) simulation2d < gridMesher
    properties (Constant)
        void_density = 1e-3;  % density of void elements
        void_density_interpolated = 1e-9;    % density of void elements for interpolation
        solid_threshold = 0.5;  % threshold for solid elements
    end
    properties(GetAccess = 'public', SetAccess = 'protected')
        m_vectorize = 0; % whether to use vectorized implementation (much faster)
        m_design;   % design variable (density)
        m_solidElems;   % solid elements
        m_solidNodes;    % solid nodes
        m_interpolation; % interpolation method
        m_materials;    % list of materials
        m_numMaterials;      % number of materials
        m_materialIndices;  % indices of materials for each element
        m_numScenarios; % number of load scenarios
    end

    methods (Abstract)
        obj = preProcess(obj)       
        obj = solve(obj)
        obj = postProcess(obj)
    end
    methods
        function obj = simulation2d(brep,numElements,materials,...
                interpolation,numScenarios,uniformGrid)
            % Constructor for the simulation2d class
            % Inputs:
            %   brep: boundary representation of the geometry
            %   numElements: number of elements in the mesh
            %   materials: list of materials for the simulation
            %   interpolation: interpolation method for design variables
            %   numScenarios: number of load scenarios
            %   uniformGrid: flag for uniform grid meshing (default: false)
            % Outputs:
            %   obj: instance of the simulation2d class

            obj = obj@gridMesher(brep,numElements,uniformGrid); % call superclass
            obj = obj.generateGrid();
            obj = obj.setDesign(obj.m_existingElems);
            obj.m_numScenarios = numScenarios;
            obj.m_interpolation = interpolation;
            obj.m_materials = materials; % list of materials
            obj.m_numMaterials = numel(materials);
            obj.m_materialIndices = ones(size(obj.m_design)); % initialize to the first material
            
        end
        %% SET DESIGN
        function obj = setDesign(obj, design)
            % Set design variable
            % Inputs:
            %   design: design variable for the simulation, can be numeric (e.g., density-based TO) or logical (e.g., level-set TO)
            % Outputs:
            %   obj: instance of the simulation2d class with updated design variable

            % Check if design has the same size as existing elements
            if ~isequal(size(design), size(obj.m_existingElems))
                error('Design must have the same size as existing elements.');
            end
            
            % Check if design values are in the range [0, 1]
            if any(design(:) < -1e-8) || any(design(:) > 1+1e-8) ...
                    || any(isnan(design(:))) || any(isinf(design(:)))
                error('Design values must be in the range [0, 1].');
            end
            % Update design variable
            obj.m_design = design;
            % Find solid elements
            obj.m_solidElems = obj.m_existingElems;
            % To avoid empty solid nodes at early optimization steps
            threshold = min(obj.solid_threshold, max(obj.m_design(:)));
            obj.m_solidElems(obj.m_design < threshold) = 0;
            % Find solid nodes
            obj.m_solidNodes = zeros(obj.m_ny+1, obj.m_nx+1);
            if (obj.m_vectorize )
                % Create logical arrays for solid elements
                solidElemsLogical = obj.m_solidElems ~= 0;
                % Add contributions from adjacent nodes
                obj.m_solidNodes(1:end-1, 1:end-1) = obj.m_solidNodes(1:end-1, 1:end-1) | solidElemsLogical;
                obj.m_solidNodes(2:end, 1:end-1) = obj.m_solidNodes(2:end, 1:end-1) | solidElemsLogical;
                obj.m_solidNodes(1:end-1, 2:end) = obj.m_solidNodes(1:end-1, 2:end) | solidElemsLogical;
                obj.m_solidNodes(2:end, 2:end) = obj.m_solidNodes(2:end, 2:end) | solidElemsLogical;
            else
                for elx = 1:obj.m_nx
                    for ely = 1:obj.m_ny
                        if (~obj.m_solidElems(ely,elx)), continue;end
                        obj.m_solidNodes(ely,elx) = 1;
                        obj.m_solidNodes(ely+1,elx) = 1;
                        obj.m_solidNodes(ely,elx+1) = 1;
                        obj.m_solidNodes(ely+1,elx+1) = 1;
                    end
                end
            end
        end

        %% COMPUTE FIELD AT NODES
        function [obj,nodalField] = computeNodalField(obj,elemField)
            % Compute nodal field from element field
            % Inputs:
            %   elemField: field values at elements (e.g., stress, strain)
            % Outputs:
            %   obj: instance of the simulation2d class with updated nodal field
            %   nodalField: field values at nodes (averaged from elements)
            % Check if elemField is a valid matrix
            if ~isnumeric(elemField) || ~ismatrix(elemField)
                error('Element field must be a numeric matrix.');
            end
            % Check if elemField has the same size as existing elements
            if ~isequal(size(elemField), size(obj.m_existingElems))
                error('Element field must have the same size as existing elements.');
            end
            % Initialize nodal field and solid neighbors
            % nodalField will store the averaged field values at nodes
            % solidNeighbors will count the number of solid elements contributing to each node
            % Initialize nodalField and solidNeighbors
            nodalField = zeros(obj.m_ny+1,obj.m_nx+1);
            solidNeighbors = zeros(obj.m_ny+1,obj.m_nx+1);
            localNodeIds = [0 0;1 0;0 1;1 1];
            for elx = 1:obj.m_nx
                for ely = 1:obj.m_ny
                    if (~obj.m_solidElems(ely,elx)),continue;end
                    val = elemField(ely,elx);
                    for i = 1:4
                        node = [ely+localNodeIds(i,1),elx+localNodeIds(i,2)];
                        solidNeighbors(node(1),node(2)) = solidNeighbors(node(1),node(2)) + 1;
                        nodalField(node(1),node(2)) = nodalField(node(1),node(2)) + val;
                    end
                end
            end
            nodalField = nodalField ./ max(1,solidNeighbors);
        end
    end
    methods(Static)
        %% GAUSSIAN QUADRATURE POINTS FOR LINE INTEGRATION
        function [xi_GQ, wt_GQ] = GaussQLine()
            % Gauss quadrature points and weights for a line element
            % Outputs:
            %   xi_GQ: Gauss quadrature points for the line element (-1 to 1)
            %   wt_GQ: Gauss quadrature weights for the line element
            % This function returns the Gauss quadrature points and weights
            % for a line element, which are used for numerical integration
            % over the element. The points are defined in the range [-1, 1]
            % and the weights are set to 1 for each point.
 
            xi_GQ = [-0.577350269189626 0.577350269189626];
            wt_GQ = [1 1];
        end
        %% LINEAR LINE SHAPE FUNCTIONS
        function [N,gradN] = edgeShapeFunction(xi)
            % Linear shape functions for a line element
            % Inputs:
            %   xi: Gauss quadrature point in the range [-1, 1]
            % Outputs:
            %   N: shape function values at the Gauss quadrature point
            %   gradN: gradient of the shape function at the Gauss quadrature point
            % This function computes the shape function values and their
            % gradients for a line element at the specified Gauss quadrature point.
            % The shape functions are linear and defined in the range [-1, 1].
            % The shape functions are defined as:
            % N1 = (1-xi)/2
            % N2 = (1+xi)/2
            % The gradients are computed as:
            % gradN1 = -1/2
            % gradN2 = 1/2

            N = [(1-xi)/2;
                (1+xi)/2];
                
            gradN = [-1/2;
                1/2];
        end
        %% GAUSSIAN QUADRATURE POINTS FOR BILINEAR QUADRILATERAL ELEMENTS
        function [xi_GQ,eta_GQ, wt_GQ] = GaussQuad()
            % Gauss quadrature points and weights for a bilinear quadrilateral element
            % Outputs:
            %   xi_GQ: Gauss quadrature points in the xi direction (-1 to 1)
            %   eta_GQ: Gauss quadrature points in the eta direction (-1 to 1)
            %   wt_GQ: Gauss quadrature weights for the quadrilateral element
            % This function returns the Gauss quadrature points and weights
            % for a bilinear quadrilateral element, which are used for numerical
            % integration over the element. The points are defined in the range
            % [-1, 1] for both xi and eta directions, and the weights are
            % set to 1 for each point.
            xi_GQ = [-1/sqrt(3) 1/sqrt(3)  1/sqrt(3) -1/sqrt(3)];
            eta_GQ = [-1/sqrt(3) -1/sqrt(3)  1/sqrt(3) 1/sqrt(3)];
            wt_GQ = [1 1 1 1];
        end
        %% QUADRILATERAL BI-LINEAR SHAPE FUNCTIONS
        function  [N,gradN] = QuadShapeFunction(xi,eta)
            % Bilinear shape functions for a quadrilateral element
            % Inputs:
            %   xi: Gauss quadrature point in the xi direction (-1 to 1)
            %   eta: Gauss quadrature point in the eta direction (-1 to 1)
            % Outputs:
            %   N: shape function values at the Gauss quadrature point
            %   gradN: gradient of the shape function at the Gauss quadrature point
            % This function computes the shape function values and their
            % gradients for a quadrilateral element at the specified Gauss
            % quadrature point. The shape functions are bilinear and defined
            % in the range [-1, 1] for both xi and eta directions.
            % The shape functions are defined as:
            % N1 = (1-xi)(1-eta)/4
            % N2 = (1+xi)(1-eta)/4
            % N3 = (1+xi)(1+eta)/4
            % N4 = (1-xi)(1+eta)/4
            % The gradients are computed as:
            % gradN1 = [-1/4*(1-eta) -1/4*(1-xi)]
            % gradN2 = [ 1/4*(1-eta) -1/4*(1+xi)]
            % gradN3 = [ 1/4*(1+eta)  1/4*(1+xi)]
            % gradN4 = [-1/4*(1+eta)  1/4*(1-xi)]
            % where xi and eta are the Gauss quadrature points in the
            % xi and eta directions, respectively.
            N = 0.25*[(1-xi)*(1-eta) (1+xi)*(1-eta) (1+xi)*(1+eta) (1-xi)*(1+eta)];
            gradN = 0.25*[eta-1 1-eta eta+1 -eta-1; xi-1 -xi-1 xi+1 1-xi];
        end
        %% JACOBIAN MATRIX
        function [J] = Jacobian(xNodes,yNodes,xi,eta)
            % Jacobian matrix for a quadrilateral element
            % Inputs:
            %   xNodes: x-coordinates of the nodes of the quadrilateral element
            %   yNodes: y-coordinates of the nodes of the quadrilateral element
            %   xi: Gauss quadrature point in the xi direction (-1 to 1)
            %   eta: Gauss quadrature point in the eta direction (-1 to 1)
            % Outputs:
            %   J: Jacobian matrix for the quadrilateral element
            % This function computes the Jacobian matrix for a quadrilateral
            % element at the specified Gauss quadrature point. The Jacobian
            % matrix is used to transform coordinates from the local element
            % space to the global coordinate system. The Jacobian is computed
            % using the shape function gradients and the coordinates of the
            % nodes of the quadrilateral element.
  
            [~,gradN] = fea2d.QuadShapeFunction(xi,eta);
            J = zeros(2,2);
            J(1,1) = gradN(1,:)*xNodes';
            J(1,2) = gradN(2,:)*xNodes';
            J(2,1) = gradN(1,:)*yNodes';
            J(2,2) = gradN(2,:)*yNodes';
        end
    end
end