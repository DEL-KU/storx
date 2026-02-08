
classdef test_parametricOptimization < matlab.unittest.TestCase
    %% Test class for parametric optimization
    methods (Test)

        function testCantileverBeamRS(testCase)
            %% Create Problem
            function fem = createProblem(brep)
                vectorize = true;
                numElements = 500; % mesh
                material.E = 100e9; material.nu = 0.3; material.rho = 1; % material
                numScenarios = 1;
                fem = fea2d_elasticity(brep,numElements,material,vectorize,numScenarios);
                fem = fem.fixEdge([2,10]);
                fem = fem.applyYForceOnEdge(6,-1e5);
            end

            %% Create Geometry from Parameters
            function geom = createGeom(params)
                % params
                L = 2.0;%length
                H = 1;%height
                h = 0.2; % length of force edge
                a = params(1); % corner cutouts
                b = params(2); % left edge cutout
                geom.vertices = [b 0;0 -b; 0 -H/2;L-a -H/2;L -H/2+a; L -h/2; L h/2; L H/2-a; L-a H/2; 0 H/2; 0 b ]';

                geom.segments = [1 1 2 0; 1 2 3 0;1 3 4 0;1 4 5 0;1 5 6 0;1 6 7 0;1 7 8 0;1 8 9 0;1 9 10 0;1 10 11 0;1 11 1 0]';
            end
            %% Problem Definition
            params0.value = [0.1 0.15];
            params0.lb = [0.05 0.05];
            params0.ub = [0.4 0.4];

            objective = 'compliance'; % objective
            constraints.area = 1.8; % constraint value
            constraints.type = 'ineq'; % constraint type: 'eq' or 'ineq'
            %% Construct Optimizer
            brepHandle = @createGeom;
            solverHandle = @createProblem;

            % Optimization method:
            %  - RS: Random Search
            %  - FD: Finite Difference
            %  - MS: Multi-Start
            %  - GS: Global Search
            method = "RS";
            exportGIF = false; % export optimization progress as GIF
            testMode = true; % test mode: no export, no plotting
            parOpt = parameterOpt2d(brepHandle,solverHandle,params0, ...
                objective,constraints, ...
                [],[],method,exportGIF,testMode);

            %% Optimize
            parOpt = parOpt.optimize();

            %% Verify convergence
            testCase.verifyTrue(ismember(parOpt.m_results.flag, [1 2 3]), ...
                 'Parametric optimization did not converge to a solution.');

        end
    end
end
