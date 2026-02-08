classdef test_fea2d_thermal < matlab.unittest.TestCase
    %% Test class for thermal FEA
    methods (Test)

        function testThermalFEA(testCase)
            %% Problem definition
            % This test checks the temperature in a vertical bar
            % under a thermal flux applied at the top edge.
            % The expected maximum temperature is 10.   
            vectorize = true;
            material.k = 1; % thermal conductivity
            numElements = 1000;
            verticalBar.vertices = [0 0; 1 0; 1 10; 0 10]';
            verticalBar.segments = [1 1 2 0 ;1 2 3 0;1 3 4 0;1 4 1 0]';
            fem = fea2d_thermal(verticalBar,numElements,material,vectorize);
            fem = fem.fixEdge(1,0);
            fem = fem.applyFlux(3,1);

            %% Solve
            fem = fem.preProcess();
            fem = fem.solve();
            fem = fem.postProcess();

            %% Check
            TMax = max(fem.m_sol(:));
            testCase.verifyEqual(TMax,10,'AbsTol',1e-2);
        end
    end
end
