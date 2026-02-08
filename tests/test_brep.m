classdef test_brep < matlab.unittest.TestCase
    %% Test class for brep2d geometries
    methods (Test)

        function testBrepArea_Rectangle(testCase)
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
                1 4 1 0]';% left edge

            % Act
            geom = brep2d(rectBrep);  % Construct brep2d object

            % Assert — B-rep vertices and segments
            testCase.verifyEqual(geom.getNumberofBrepVertices(), 4, ...
                'Should have 4 vertices.');
            testCase.verifyEqual(geom.getNumberofBrepSegments(), 4, ...
                'Should have 4 segments.');
            testCase.verifyEqual(geom.getBrepVertices(), rectBrep.vertices, ...
                'Vertices mismatch.');

            % Assert — Bounding box
            testCase.verifyEqual(geom.getBoundingBox(), [0 2;0 1], ...
                'Bounding box mismatch.');
            testCase.verifyEqual(geom.getBoxSizes(), [2 1 1], ...
                'Box sizes mismatch.');
            testCase.verifyEqual(round(geom.getCenter(),3), [1;0.5], ...
                'Center mismatch.');

            % Assert — Area & perimeter
            expectedArea = 2*1;
            expectedPerimeter = 2*(2+1);
            testCase.verifyEqual(round(geom.getArea(),3), expectedArea, ...
                'Area mismatch.');
            testCase.verifyEqual(round(geom.getPerimeter(),3), expectedPerimeter, ...
                'Perimeter mismatch.');

            % Assert — BRepArea calculation matches getArea()
            testCase.verifyEqual(round(geom.brepArea(),3), expectedArea, ...
                'BRep area mismatch.');
        end


        function testBeamWithChamfersLoads(testCase)
            % Arrange
            filename = 'BeamWithChamfers.brep';

            vertices = [ ...
                0    1.7  2    2    1.7  0;    % x-coordinates
                0    0    0.3  0.7  1    1];   % y-coordinates

            % Act
            geom = brep2d(filename);  % Load from file

            % Assert — B-rep vertices and segments
            testCase.verifyEqual(geom.getNumberofBrepVertices(), 6, ...
                'Should have 6 vertices.');
            testCase.verifyEqual(geom.getNumberofBrepSegments(), 6, ...
                'Should have 6 segments.');
            testCase.verifyEqual(geom.getBrepVertices(), vertices, ...
                'Vertices mismatch.');

            % Assert — Bounding box
            testCase.verifyEqual(geom.getBoundingBox(), [0 2;0 1], ...
                'Bounding box mismatch.');
            testCase.verifyEqual(geom.getBoxSizes(), [2 1 1], ...
                'Box sizes mismatch.');

            expectedCenter = mean(vertices, 2);
            testCase.verifyEqual(round(geom.getCenter(),3), round(expectedCenter,3), ...
                'Center mismatch.');

            % Assert — Area & perimeter
            expectedArea = round(2-0.3^2,2);
            expectedPerimeter = round(1.7+1.7+0.4+1+2*(0.3*2/sqrt(2)),2);
            testCase.verifyEqual(round(geom.getArea(),2), expectedArea, ...
                'Area mismatch.');
            testCase.verifyEqual(round(geom.getPerimeter(),2), expectedPerimeter, ...
                'Perimeter mismatch.');
        end

        function testBeamWithFilletsConstructed(testCase)
            % Arrange
            BeamWithFillets.vertices = [0 0; 1.7 0; 2 0.3; 2 0.7; 1.7 1; 0 1; 1.7 0.3; 1.7 0.7]';
            BeamWithFillets.segments = [ ...
                1 1 2 0; ...
                2 2 3 -7; ...
                1 3 4 0; ...
                2 4 5 -8; ...
                1 5 6 0; ...
                1 6 1 0 ]';

            % Act
            geom = brep2d(BeamWithFillets);

            % Assert
            testCase.verifyEqual(geom.getNumberofBrepVertices(), 8, ...
                'Should have 8 vertices.');
            testCase.verifyEqual(geom.getNumberofBrepSegments(), 6, ...
                'Should have 6 segments.');


            expectedArea = round(1.7+0.4*0.3+pi*0.3^2/2,2);
            expectedPerimeter = round(1.7+1.7+0.4+1+(pi*0.3),2);

            testCase.verifyEqual(round(geom.getArea(), 2), expectedArea, ...
                'Area mismatch.');
            testCase.verifyEqual(round(geom.getPerimeter(), 2), expectedPerimeter, ...
                'Perimeter mismatch.');

            expectedBox = [0 2; 0 1];
            testCase.verifyEqual(round(geom.getBoundingBox(),3), expectedBox, ...
                'Bounding box mismatch.');

            testCase.verifyEqual(geom.getBoxSizes(), [2 1 1], ...
                'Box size mismatch.');

            expectedCenter = mean(BeamWithFillets.vertices,2);
            testCase.verifyEqual(round(geom.getCenter(),3), round(expectedCenter,3), ...
                'Center mismatch.');

            testCase.verifyEqual(size(geom.getBrepVertices(),2), 8, ...
                'Vertex count mismatch.');

            testCase.verifyEqual(length(geom.getSegLengths()), 6, ...
                'Segment lengths count mismatch.');
        end

        function testBeamWithHoleConstructed(testCase)
            % Arrange
            BeamWithHole.vertices = [0 0;1 0;1 0.3;1 0.7;2 0;2 1;0 1;1 0.5]';
            BeamWithHole.segments = [1 1 2 0; -1 2 3 0; 2 3 4 8; 2 4 3 8; -1 3 2 0; 1 2 5 0; 1 5 6 0; 1 6 7 0; 1 7 1 0]';

            % Act
            geom = brep2d(BeamWithHole);

            % Assert
            testCase.verifyEqual(geom.getNumberofBrepVertices(), 8, ...
                'Should have 8 vertices.');
            testCase.verifyEqual(geom.getNumberofBrepSegments(), 9, ...
                'Should have 9 segments.');


            expectedArea = round(2-pi*0.2^2,2);
            expectedPerimeter = round(6+2*(pi*0.2),2);

            testCase.verifyEqual(round(geom.getArea(), 2), expectedArea, ...
                'Area mismatch.');
            testCase.verifyEqual(round(geom.getPerimeter(), 2), expectedPerimeter, ...
                'Perimeter mismatch.');

            expectedBox = [0 2; 0 1];
            testCase.verifyEqual(round(geom.getBoundingBox(),3), expectedBox, ...
                'Bounding box mismatch.');

            testCase.verifyEqual(geom.getBoxSizes(), [2 1 1], ...
                'Box size mismatch.');

            expectedCenter = mean(BeamWithHole.vertices,2);
            testCase.verifyEqual(round(geom.getCenter(),3), round(expectedCenter,3), ...
                'Center mismatch.');

            testCase.verifyEqual(size(geom.getBrepVertices(),2), 8, ...
                'Vertex count mismatch.');

            testCase.verifyEqual(length(geom.getSegLengths()), 9, ...
                'Segment lengths count mismatch.');
        end

        function testAnnulusConstructed(testCase)
            % Act
            filename = 'Annulus.brep';

            vertices = [ ...
                -1    0    1    0.5   0   -0.5   0    0    0.0 ;  % x-coordinates
                0   -1    0    0     -0.5  0    0.5  1    0.0 ]; % y-coordinates


            geom = brep2d(filename);

            % Assert
            testCase.verifyEqual(geom.getNumberofBrepVertices(), 9, ...
                'Should have 9 vertices.');
            testCase.verifyEqual(geom.getNumberofBrepSegments(), 10, ...
                'Should have 10 segments.');


            expectedArea = round(pi*(1-0.5^2),2);
            expectedPerimeter = round(2*pi*(1+0.5),2);


            testCase.verifyLessThan(abs(geom.getArea()-expectedArea)/expectedArea, 0.01, ...
                'Area mismatch.');


            testCase.verifyEqual(round(geom.getPerimeter(), 2), expectedPerimeter, ...
                'Perimeter mismatch.');

            expectedBox = [-1 1; -1 1];
            testCase.verifyEqual(round(geom.getBoundingBox(),3), expectedBox, ...
                'Bounding box mismatch.');

            testCase.verifyEqual(geom.getBoxSizes(), [2 2 1], ...
                'Box size mismatch.');

            expectedCenter = mean(vertices,2);
            testCase.verifyEqual(round(geom.getCenter(),3), round(expectedCenter,3), ...
                'Center mismatch.');

            testCase.verifyEqual(size(geom.getBrepVertices(),2), size(vertices,2), ...
                'Vertex count mismatch.');

            testCase.verifyEqual(length(geom.getSegLengths()), 10, ...
                'Segment lengths count mismatch.');
        end


    end
end
