classdef test_topopt2d_levelset < matlab.unittest.TestCase
    %% Test class for levelset2d topopt
    methods (Test)

        function testLbrackt_levelset_SO_without_holes(testCase)
            %% General Parameters
            vectorize = true;
            uniformGrid = 1; % needed for the Hamilton-Jacobi solver
            exportGIF = false;

            %% Solvers
            feaClass = @fea2d_elasticity;
            shapeoptClass = @standardHJ2d_elasticity;

            %% Optimizer Parameters
            interpolation = 'none';
            maxNumIters = 20;
            penaltyStruct = struct('min',1,'max',1,'inc',0);

            %% Problem definition
            brep = 'LBracketNoFillet.brep'; % geometry
            numElements = 3200; % mesh
            material.E = 100e9; material.nu = 0.3; material.rho = 1; % material
            numScenarios = 1;

            %% Construct FEA Solver
            solver = feaClass(brep,numElements,material,vectorize,numScenarios, ...
                interpolation,penaltyStruct,uniformGrid); % call superclass

            solver = solver.fixEdge(6);
            solver = solver.applyYForceOnEdge(3,-1e5);

            solver = solver.preProcess(); % FEA pre-processing

            %% Objective and Constraints
            objective = standardHJComplianceElasticity(solver);

            volumeFraction = 0.5;
            constraints = {volume(solver, volumeFraction)};

            % manufacturing constraints
            mfgConstraints = {minimumFeatureSize_conv(solver)
                retain_levelset(solver,3)};

            testMode = true; % test mode: no export, no plotting
            %% Construct Optimizer
            nHolesX = 0; nHolesY = 0; r0 = 0;

            shapeopt = shapeoptClass(solver, ...
                objective,constraints,mfgConstraints, ...
                nHolesX,nHolesY,r0, ...
                maxNumIters,exportGIF,testMode);

            %% Optimize
            shapeopt = shapeopt.optimize();

            %% Verify convergence
            testCase.verifyTrue(ismember(shapeopt.m_flag, [1 2]), ...
                'LevelSet Shape Optimization without hole initialization for elasticity did not converge to a solution.');

        end


        function testLbrackt_levelset_SO_with_holes(testCase)
            %% General Parameters
            vectorize = true;
            uniformGrid = 1; % needed for the Hamilton-Jacobi solver
            exportGIF = false;

            %% Solvers
            feaClass = @fea2d_elasticity;
            shapeoptClass = @standardHJ2d_elasticity;

            %% Optimizer Parameters
            interpolation = 'none';
            maxNumIters = 20;
            penaltyStruct = struct('min',1,'max',1,'inc',0);

            %% Problem definition
            brep = 'LBracketNoFillet.brep'; % geometry
            numElements = 3200; % mesh
            material.E = 100e9; material.nu = 0.3; material.rho = 1; % material
            numScenarios = 1;

            %% Construct FEA Solver
            solver = feaClass(brep,numElements,material,vectorize,numScenarios, ...
                interpolation,penaltyStruct,uniformGrid); % call superclass

            solver = solver.fixEdge(6);
            solver = solver.applyYForceOnEdge(3,-1e5);

            solver = solver.preProcess(); % FEA pre-processing

            %% Objective and Constraints
            objective = standardHJComplianceElasticity(solver);

            volumeFraction = 0.5;
            constraints = {volume(solver, volumeFraction)};

            % manufacturing constraints
            mfgConstraints = {minimumFeatureSize_conv(solver)
                retain_levelset(solver,3)};

            testMode = true; % test mode: no export, no plotting
            %% Construct Optimizer
            nHolesX = 3; nHolesY = 3; r0 = 0.5;

            shapeopt = shapeoptClass(solver, ...
                objective,constraints,mfgConstraints, ...
                nHolesX,nHolesY,r0, ...
                maxNumIters,exportGIF,testMode);

            %% Optimize
            shapeopt = shapeopt.optimize();

            %% Verify convergence
            testCase.verifyTrue(ismember(shapeopt.m_flag, [1 2]), ...
                'LevelSet Shape Optimization with hole initialization for elasticity did not converge to a solution.');

        end

        function testLbrackt_levelset_TO(testCase)
            %% General Parameters
            vectorize = true;
            uniformGrid = 1; % needed for the Hamilton-Jacobi solver
            exportGIF = false;

            %% Solvers
            feaClass = @fea2d_elasticity;
            topoptClass = @modifiedHJ2d_elasticity;

            %% Optimizer Parameters
            interpolation = 'none';
            maxNumIters = 20;
            penaltyStruct = struct('min',1,'max',1,'inc',0);

            %% Problem definition
            brep = 'LBracketNoFillet.brep'; % geometry
            numElements = 3200; % mesh
            material.E = 100e9; material.nu = 0.3; material.rho = 1; % material
            numScenarios = 1;

            %% Construct FEA Solver
            solver = feaClass(brep,numElements,material,vectorize,numScenarios, ...
                interpolation,penaltyStruct,uniformGrid); % call superclass

            solver = solver.fixEdge(6);
            solver = solver.applyYForceOnEdge(3,-1e5);

            solver = solver.preProcess(); % FEA pre-processing

            %% Objective and Constraints
            objective = modifiedHJComplianceElasticity(solver);

            volumeFraction = 0.5;
            constraints = {volume(solver, volumeFraction)};

            % manufacturing constraints
            mfgConstraints = {minimumFeatureSize_conv(solver)
                retain_levelset(solver,3) };

            testMode = true; % test mode: no export, no plotting
            %% Construct Optimizer
            topWeight = 10;
            topopt = topoptClass(solver, ...
                objective,constraints,mfgConstraints, ...
                topWeight, ...
                maxNumIters,exportGIF,testMode);

            %% Optimize
            topopt = topopt.optimize();

            %% Verify convergence
            testCase.verifyTrue(ismember(topopt.m_flag, [1 2]), ...
                'LevelSet Topology Optimization for elasticity did not converge to a solution.');

        end

    end
end
