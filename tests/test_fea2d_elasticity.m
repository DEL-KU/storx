classdef test_fea2d_elasticity < matlab.unittest.TestCase
    %% Test class for 2D elasticity FEA
    methods (Test)

        function testElasticityFEA_LoadIntegration(testCase)
            %% Test for 2D elasticity FEA with load integration
            % This test checks the nodal force distribution on a loaded edge
            % of a cantilever beam under a unit traction force applied on the
            % right edge.
            % The total force on the edge should equal the applied force.
            %% Parameters
            vectorize = true;
            material.E = 2e11; material.nu = 0.28; material.rho = 1;
            numElements = 1000;
            H = 1;
            cantileverBeam.vertices = [0 0; 2 0; 2 H; 0 H]';
            cantileverBeam.segments = [1 1 2 0 ;1 2 3 0;1 3 4 0;1 4 1 0]';

            fem = fea2d_elasticity(cantileverBeam,numElements,material,vectorize);
            fem = fem.fixEdge(4);
            force = 1;           % left edge fixed
            fem = fem.applyXForceOnEdge(3,force);   % unit traction on right edge

            fem = fem.assembleBC();
            f = full(fem.m_f);

            %% Total force check
            totalFx = sum(f(1:2:end));
            testCase.verifyEqual(totalFx, force, 'AbsTol',1e-12)

            %% Nodal force distribution on loaded edge
            loadedSegments = find(fem.m_edges(5,:) == 3);   % segments belonging to edge #3
            allNodes = unique(fem.m_edges(1:fem.m_nodesPerEdge, loadedSegments)); % node numbers
            xDOF = 2*allNodes-1;
            fxNodes = f(xDOF);

            %% Check that all nodal x-forces sum to totalFx
            testCase.verifyEqual(sum(fxNodes), totalFx, 'AbsTol',1e-12)

            if fem.m_nodesPerEdge == 2
                numNodesOnEdge = numel(fxNodes);
                numSegmentsOnEdge = numNodesOnEdge-1;
                segForce = 1 * (1/numSegmentsOnEdge); % per-segment force

                expectedFxNodes = zeros(numNodesOnEdge,1);
                expectedFxNodes(1)   = segForce/2;
                expectedFxNodes(end) = segForce/2;
                expectedFxNodes(2:end-1) = segForce;  % interior nodes
                testCase.verifyEqual(fxNodes, expectedFxNodes, ...
                    'AbsTol',1e-12, ...
                    'Nodal x-forces not distributed as expected')
            end

        end


        function testElasticityFEA_Tension(testCase)
            %% Test for 2D elasticity FEA with tension in y-direction
            % This test checks the displacement and stress in a vertical bar
            % under a vertical force applied at the top edge.
            %% General parameters
            vectorize = true;

            %% Problem definition
            material.E = 2e11;  material.nu = 0.28;
            material.rho = 1;
            numElements = 1000;
            verticalBar.vertices = [0 0; 1 0; 1 10; 0 10]';
            verticalBar.segments = [1 1 2 0 ;1 2 3 0;1 3 4 0;1 4 1 0]';
            fem = fea2d_elasticity(verticalBar,numElements,material,vectorize);
            fem = fem.fixYOfEdge(1);
            fem = fem.applyYForceOnEdge(3,100);
            %% Solve
            fem = fem.preProcess();
            fem = fem.solve();
            fem = fem.postProcess(true);
            %% Test assertions
            % Check displacement in y-direction
            expectedDisp = 5e-9;
            diffDisp = abs(max(fem.m_sol(2:2:end)) - expectedDisp);  % absolute diff
            testCase.verifyLessThan(diffDisp, 1e-12, ...
                'Max displacement in y-direction is not within 1e-12 of the expected value.');

            % Assume you want the y-direction normal stress = sigma_yy = (2,2) component
            scenarioId = 1;  % 1st scenario
            expectedStress = 100;   % Pa

            % Extract sigma_yy for all elements in this scenario
            sigma_yy = fem.m_stressTensor(:,2,2,scenarioId);

            % Compute absolute difference
            diffStress = abs(max(sigma_yy) - expectedStress);

            % Verify the max absolute difference is within tolerance
            testCase.verifyLessThan(diffStress, 1e-9, ...
                'Max sigma_yy is not within 1e-9 of the expected value.');

            % Check sigma_xx across all elements is approximately zero
            sigma_xx = fem.m_stressTensor(:,1,1,scenarioId);
            diffXX = max(abs(sigma_xx));
            testCase.verifyLessThan(diffXX, 1e-9, ...
                'Stress in x-direction is not within 1e-9 of zero for tension in y-direction.');

            % Check tau_xy (shear) across all elements is approximately zero
            tau_xy = fem.m_stressTensor(:,1,2,scenarioId);
            diffXY = max(abs(tau_xy));
            testCase.verifyLessThan(diffXY, 1e-9, ...
                'Shear stress is not within 1e-9 of zero for tension in y-direction.');
        end

        function testElasticityFEA_Bending(testCase)
            %% Test for 2D elasticity FEA with bending in y-direction
            % This test checks the displacement and stress in a vertical bar
            % under a vertical force applied at the top edge.
            %% General parameters
            vectorize = true;

            %% Problem definition
            material.E = 2e11;  material.nu = 0.; material.rho = 1;
            numElements = 10000;
            verticalBar.vertices = [0 0; 1 0; 1 10; 0 10]';
            verticalBar.segments = [1 1 2 0 ;1 2 3 0;1 3 4 0;1 4 1 0]';
            fem = fea2d_elasticity(verticalBar,numElements,material,vectorize);
            fem = fem.fixEdge(1);
            fem = fem.applyXForceOnEdge(3,100);
            %% Solve
            fem = fem.preProcess();
            fem = fem.solve();
            fem = fem.postProcess(true);
            %% Test assertions
            % Check displacement in y-direction
            expectedDef = 2e-6;    % m

            diffDef = abs(max(fem.m_def(:)) - expectedDef);  % absolute diff
            testCase.verifyLessThan(diffDef, 1e-5, ...
                'Max displacement in y-direction is not within 1e-8 of the expected value.');

            % Assume you want the x-direction normal stress = sigma_xx = (1,1) component
            scenarioId = 1;  % 1st scenario
            expectedStress = 6000;   % Pa

            maxStress = max(fem.m_vonMisesElems(:,scenarioId));
            diffStress = max(abs(maxStress - expectedStress));
            testCase.verifyLessThan(diffStress, 100, ...
                'Max von mises stress is not within 100 Pa of the expected value.');
        end

    end
end
