classdef test_fea2d_thermoelasticity < matlab.unittest.TestCase
    %% Test class for thermoelastic FEA
    methods (Test)

        function testThermoelasticityFEA(testCase)
            %% Problem definition
            % This test checks the temperature and displacement in a vertical bar
            % under a thermal flux applied at the top edge.
            % The expected maximum temperature is 66.667 (C) and the maximum deformation is 0.001 (m).
            
            vectorize = true;
            %% Problem definition
            material.E = 2e11;  material.nu = 0; material.rho = 7850;
            material.k = 45;
            material.alpha = 1e-5;
            TReference = 0;

            numElements = 100;
            xmax = 3;
            domain.vertices = [0 0; xmax 0; xmax 1; 0 1]';
            domain.segments = [1 1 2 0 ;1 2 3 0;1 3 4 0;1 4 1 0]';

            fem_thermal = fea2d_thermal(domain,numElements,material,vectorize);
            fem_thermal = fem_thermal.fixEdge(4,0); % fixed T
            Q = 1000;
            fem_thermal = fem_thermal.applyFlux(2,Q); % flux

            tau = 0;
            fem_elasticity = fea2d_elasticity(domain,numElements,material,vectorize);
            fem_elasticity = fem_elasticity.fixXOfEdge(4);
            fem_elasticity = fem_elasticity.fixYOfEdge([1,3]);
            fem = fea2d_thermoelasticity(fem_thermal,fem_elasticity, ...
                TReference,domain,numElements,material,vectorize);

            % Analytical solution:
            E = material.E; nu = material.nu; k = material.k;alpha = material.alpha;
            lambda = (E*nu)/((1+nu)*(1-2*nu));
            mu = E/(2*(1+nu));

            u_max = Q/k*alpha * (lambda+mu)/(lambda+2*mu) * xmax^2 + tau/(lambda+2*mu)*xmax;
            T_max = Q/k*xmax;

            %% Solve
            fem = fem.preProcess();
            fem = fem.solve();
            fem = fem.postProcess();

            %% Check
            uMax = max(fem.m_elasticitySolver.m_def(:));
            TMax = max(fem.m_thermalSolver.m_sol(:));
            testCase.verifyEqual(uMax,u_max,'AbsTol',1e-2);
            testCase.verifyEqual(TMax,T_max,'AbsTol',1e-2);
        end
    end
end
