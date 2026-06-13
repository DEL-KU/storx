
classdef test_parametricOptimization < matlab.unittest.TestCase
    %% Test class for parametric optimization

    % Shared geometry and problem setup helpers
    methods (Static, Access = private)

        function fem = createProblem(brep)
            vectorize = true;
            numElements = 500;
            material.E = 100e9; material.nu = 0.3; material.rho = 1;
            numScenarios = 1;
            fem = fea2d_elasticity(brep, numElements, material, vectorize, numScenarios);
            fem = fem.fixEdge([2, 10]);
            fem = fem.applyYForceOnEdge(6, -1e5);
        end

        function geom = createGeom(params)
            L = 2.0; H = 1; h = 0.2;
            a = params(1);
            b = params(2);
            geom.vertices = [b 0; 0 -b; 0 -H/2; L-a -H/2; L -H/2+a; ...
                             L -h/2; L h/2; L H/2-a; L-a H/2; 0 H/2; 0 b]';
            geom.segments = [1 1 2 0; 1 2 3 0; 1 3 4 0; 1 4 5 0; 1 5 6 0; ...
                             1 6 7 0; 1 7 8 0; 1 8 9 0; 1 9 10 0; 1 10 11 0; 1 11 1 0]';
        end

        function [params0, objective, constraints] = commonProblem()
            params0.value = [0.1 0.15];
            params0.lb    = [0.05 0.05];
            params0.ub    = [0.4  0.4 ];
            objective          = 'compliance';
            constraints.area   = 1.8;
            constraints.type   = 'ineq';
        end

    end

    methods (Test)

        function testCantileverBeamRS(testCase)
            [params0, objective, constraints] = ...
                test_parametricOptimization.commonProblem();

            brepHandle   = @test_parametricOptimization.createGeom;
            solverHandle = @test_parametricOptimization.createProblem;

            exportGIF = false;
            testMode  = true;

            parOpt = parameterOpt2d_RS(brepHandle, solverHandle, params0, ...
                objective, constraints, exportGIF, testMode);

            parOpt = parOpt.optimize();

            testCase.verifyTrue(parOpt.m_results.flag >= 1, ...
                'RS: parametric optimization did not converge to a feasible solution.');
        end

        function testCantileverBeamFD(testCase)
            [params0, objective, constraints] = ...
                test_parametricOptimization.commonProblem();

            brepHandle   = @test_parametricOptimization.createGeom;
            solverHandle = @test_parametricOptimization.createProblem;

            terminationTolerance     = 1e-4;
            finiteDifferenceStepSize = 1e-4;
            exportGIF = false;
            testMode  = true;

            parOpt = parameterOpt2d_FD(brepHandle, solverHandle, params0, ...
                objective, constraints, ...
                terminationTolerance, finiteDifferenceStepSize, ...
                exportGIF, testMode);

            parOpt = parOpt.optimize();

            testCase.verifyTrue(ismember(parOpt.m_results.flag, [1 2 3]), ...
                'FD: parametric optimization did not converge.');
        end

        function testCantileverBeamMS(testCase)
            [params0, objective, constraints] = ...
                test_parametricOptimization.commonProblem();

            brepHandle   = @test_parametricOptimization.createGeom;
            solverHandle = @test_parametricOptimization.createProblem;

            terminationTolerance     = 1e-4;
            finiteDifferenceStepSize = 1e-4;
            exportGIF = false;
            testMode  = true;

            parOpt = parameterOpt2d_MS(brepHandle, solverHandle, params0, ...
                objective, constraints, ...
                terminationTolerance, finiteDifferenceStepSize, ...
                exportGIF, testMode);

            % Reduce local starts to 2 to keep the test fast
            parOpt = parOpt.setNumberOfMultiStartLocalProblems(2);

            parOpt = parOpt.optimize();

            testCase.verifyTrue(parOpt.m_results.flag >= 1, ...
                'MS: parametric optimization did not converge.');
        end

        function testCantileverBeamGS(testCase)
            [params0, objective, constraints] = ...
                test_parametricOptimization.commonProblem();

            brepHandle   = @test_parametricOptimization.createGeom;
            solverHandle = @test_parametricOptimization.createProblem;

            terminationTolerance     = 1e-4;
            finiteDifferenceStepSize = 1e-4;
            exportGIF = false;
            testMode  = true;

            parOpt = parameterOpt2d_GS(brepHandle, solverHandle, params0, ...
                objective, constraints, ...
                terminationTolerance, finiteDifferenceStepSize, ...
                exportGIF, testMode);

            parOpt = parOpt.optimize();

            testCase.verifyTrue(parOpt.m_results.flag >= 1, ...
                'GS: parametric optimization did not converge.');
        end

    end
end
