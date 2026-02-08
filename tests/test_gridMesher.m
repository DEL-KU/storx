classdef test_gridMesher < matlab.unittest.TestCase
    %% Test class for grid mesh generation
    methods (Test)

        function testGridMesh_Rectangle(testCase)
            % Arrange
            % Define a 2 x N array of vertices for a 2x1 rectangle
            rectBrep.vertices = [ ...
                0 2 2 0;  % x
                0 0 1 1]; % y

            % Segments: straight lines between successive vertices
            % segs(i,:) = [type, vs, ve, ctrl]; for straight segments ctrl=0
            rectBrep.segments = [ ...
                1 1 2 0;  % bottom edge
                1 2 3 0;  % right edge
                1 3 4 0;  % top edge
                1 4 1 0]'; % left edge

            % Act
            numElements = 41*21;   % Total number of elements
            uniformGrid = 0;
            mesh = gridMesher(rectBrep, numElements, uniformGrid);  % Construct grid
            % Generate the grid mesh
            mesh = mesh.generateGrid();  % Generate the grid mesh
            % Assert basic properties
            testCase.verifyEqual(mesh.m_nx, 41, 'm_nx should be 41.');
            testCase.verifyEqual(mesh.m_ny, 21, 'm_ny should be 21.');
            testCase.verifyEqual(mesh.m_numElems, mesh.m_nx * mesh.m_ny, ...
                'Total elements should match nx * ny.');
            testCase.verifyEqual(mesh.m_numNodes, (mesh.m_nx + 1)*(mesh.m_ny + 1), ...
                'Total nodes should match (nx+1)*(ny+1).');

            % Verify uniform grid spacing
            testCase.verifyEqual(round(mesh.m_hx,3), round(2 / mesh.m_nx,3), 'hx should match width/numX.');
            testCase.verifyEqual(round(mesh.m_hy,3), round(1 / mesh.m_ny,3), 'hy should match height/numY.');
            testCase.verifyEqual(mesh.m_ve, mesh.m_hx * mesh.m_hy, ...
                'Element area ve should be hx * hy.');

            % Check uniformGrid flag
            testCase.verifyEqual(mesh.m_uniformGrid, uniformGrid, ...
                'UniformGrid flag should match input.');

            % Check if nodes are on the boundary
            testCase.verifyTrue(any(mesh.m_isNodeOnBoundary), ...
                'There should be boundary nodes.');

            % Check sizes of coordinate arrays
            testCase.verifySize(mesh.m_nodeCoords, [2, mesh.m_numNodes], ...
                'Node coords must have 2 rows and numNodes columns.');
            testCase.verifySize(mesh.m_elemCoords, [2,mesh.m_numElems], ...
                'Elem coords must have numElems rows and 2 columns.');

            % Check some expected properties
            testCase.verifyGreaterThan(mesh.m_numExistingElems, 0, ...
                'Number of existing elements must be positive.');

            % Verify consistency of other fields
            testCase.verifyEqual(numel(mesh.m_existingElems), mesh.m_numElems, ...
                'ExistingElems should match numElems.');
            testCase.verifyEqual(numel(mesh.m_existingNodes), mesh.m_numNodes, ...
                'ExistingNodes should match numNodes.');

            % Check if edges and quad connectivities have reasonable sizes
            testCase.verifyEqual(size(mesh.m_q,2), mesh.m_numElems, ...
                'Quad element list must have one column per element.');
            testCase.verifyTrue(~isempty(mesh.m_edges), ...
                'Edges must not be empty.');
        end

        function testGridMesh_UniformGrid(testCase)
            rectBrep.vertices = [0 2 2 0; 0 0 1 1];
            rectBrep.segments = [1 1 2 0; 1 2 3 0; 1 3 4 0; 1 4 1 0]';
            numElements = 40; uniformGrid = 1;
            mesh = gridMesher(rectBrep, numElements, uniformGrid).generateGrid();

            % Check uniform spacing
            diffValue_x = abs(diff(unique(mesh.m_nodeCoords(1,:))));
            testCase.verifyLessThan(std(diffValue_x), 1e-9, 'Difference must be smaller than 1e-9.');

            diffValue_y = abs(diff(unique(mesh.m_nodeCoords(2,:))));
            testCase.verifyLessThan(std(diffValue_y), 1e-9, 'Difference must be smaller than 1e-9.');
        end

        function testGridMesh_LBracket(testCase)
            % Act
            numElements = 900;   % Total number of elements
            uniformGrid = 0;
            mesh = gridMesher('LBracketNoFillet.brep', numElements, uniformGrid);  % Construct grid
            % Generate the grid mesh
            mesh = mesh.generateGrid();  % Generate the grid mesh
            % Assert basic properties
            testCase.verifyEqual(mesh.m_numElems, mesh.m_nx * mesh.m_ny, ...
                'Total elements should match nx * ny.');
            testCase.verifyEqual(mesh.m_numNodes, (mesh.m_nx + 1)*(mesh.m_ny + 1), ...
                'Total nodes should match (nx+1)*(ny+1).');

            % Verify uniform grid spacing
            testCase.verifyEqual(round(mesh.m_hx,3), round(1 / mesh.m_nx,3), 'hx should match width/numX.');
            testCase.verifyEqual(round(mesh.m_hy,3), round(1 / mesh.m_ny,3), 'hy should match height/numY.');
            testCase.verifyEqual(mesh.m_ve, mesh.m_hx * mesh.m_hy, ...
                'Element area ve should be hx * hy.');

            % Check uniformGrid flag
            testCase.verifyEqual(mesh.m_uniformGrid, uniformGrid, ...
                'UniformGrid flag should match input.');

            % Check if nodes are on the boundary
            testCase.verifyTrue(any(mesh.m_isNodeOnBoundary), ...
                'There should be boundary nodes.');

            % Check sizes of coordinate arrays
            testCase.verifySize(mesh.m_nodeCoords, [2, mesh.m_numNodes], ...
                'Node coords must have 2 rows and numNodes columns.');
            testCase.verifySize(mesh.m_elemCoords, [2,mesh.m_numElems], ...
                'Elem coords must have numElems rows and 2 columns.');

            % Check some expected properties
            testCase.verifyGreaterThan(mesh.m_numExistingElems, 0, ...
                'Number of existing elements must be positive.');

            % Verify consistency of other fields
            testCase.verifyEqual(numel(mesh.m_existingElems), mesh.m_numElems, ...
                'ExistingElems should match numElems.');
            testCase.verifyEqual(numel(mesh.m_existingNodes), mesh.m_numNodes, ...
                'ExistingNodes should match numNodes.');

            % Check if edges and quad connectivities have reasonable sizes
            testCase.verifyEqual(size(mesh.m_q,2), mesh.m_numExistingElems, ...
                'Quad element list must have one column per element.');
            testCase.verifyTrue(~isempty(mesh.m_edges), ...
                'Edges must not be empty.');
        end
    end
end
