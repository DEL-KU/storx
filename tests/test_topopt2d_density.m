
classdef test_topopt2d_density < matlab.unittest.TestCase
    %% Test class for density topopt
    methods (Test)

        function test_elasticity_Lbracket(testCase)
            %% Solvers
            feaClass = @fea2d_elasticity;
            topoptClass = @density2d_elasticity;
            %% General Parameters
            vectorize = true;
            %% Optimizer Parameters
            interpolation = 'simp';
            update = 'OC';
            maxNumIters = 300;
            penaltyStruct = struct('min',3,'max',3,'inc',0);
            exportGIF = false; % export optimization progress as GIF
            %% Problem definition
            brep = 'LBracketNoFillet.brep'; % geometry
            numElements = 3200; % mesh
            material.E = 100e9; material.nu = 0.3; material.rho = 1; % material
            numScenarios = 1;

            %% Construct FEA Solver
            solver = feaClass(brep,numElements,material,vectorize,numScenarios, ...
                interpolation,penaltyStruct); % call superclass

            solver = solver.fixEdge(6);
            solver = solver.applyYForceOnEdge(3,-1e5);

            solver = solver.preProcess(); % FEA pre-processing

            %% Objective and Constraints
            objective = densityComplianceElasticity(solver);

            volumeFraction = 0.5;
            constraints = {volume(solver, volumeFraction)};

            % manufacturing constraints
            rmin = 1.5;
            mfgConstraints = {
                minimumFeatureSize_dist(solver, rmin)
                };

            %% Construct Optimizer
            testMode = true; % test mode: no export, no plotting
            topopt = topoptClass(solver, ...
                objective,constraints,mfgConstraints, ...
                update, ...
                maxNumIters,exportGIF, testMode);

            %% Optimize
            topopt = topopt.optimize();

            %% Verify convergence
            testCase.verifyTrue(ismember(topopt.m_flag, [1 2]), ...
                'Density TopOpt for elasticity did not converge to a solution.');
        end

        function test_thermal_squareFlux(testCase)
            %% Solvers
            feaClass = @fea2d_thermal;
            topoptClass = @density2d_thermal;

            %% General Parameters
            vectorize = true;
            exportGIF = false;
            %% Optimizer Parameters
            interpolation = 'simp';
            update = 'OC';
            maxNumIters = 500;
            penaltyStruct = struct('min',3,'max',3,'inc',0);

            %% Problem Definition
            brep = 'SquareSplit.brep'; % geometry
            numElements = 10000; % mesh
            numScenarios = 1; % # loading scenarios
            material.k = 1;
            %% Construct FEA Solver
            solver = feaClass(brep,numElements,material,vectorize,numScenarios, ...
                interpolation,penaltyStruct); % call superclass

            solver = solver.fixEdge(2,0);
            solver = solver.applyFlux(5,1);

            solver = solver.preProcess(); % FEA pre-processing

            %% Objective and Constraints
            objective = densityComplianceThermal(solver);

            volumeFraction = 0.5;
            constraints = {volume(solver, volumeFraction)};

            % manufacturing constraints
            rmin = 1.5;
            mfgConstraints = {
                minimumFeatureSize_dist(solver, rmin)
                };

            %% Construct Optimizer
            testMode = true; % test mode: no export, no plotting
            topopt = topoptClass(solver, ...
                objective,constraints,mfgConstraints, ...
                update, ...
                maxNumIters,exportGIF,testMode);

            %% Optimize
            topopt = topopt.optimize();
            %% Verify convergence
            testCase.verifyTrue(ismember(topopt.m_flag, [1 2]), ...
                'Density TopOpt for thermal did not converge to a solution.');
        end

        function test_fluid_pipeBend(testCase)
            %% Solvers
            feaClass = @fea2d_fluid;
            topoptClass = @density2d_fluid;

            %% General Parameters
            exportGIF = false;

            %% Optimizer Parameters
            interpolation = 'simp';
            update = 'MMA';
            maxNumIters = 20;

            %% Problem Definition
            brep = 'PipeBend.brep'; % geometry
            numElements = 500; % mesh
            numScenarios = 1; % # loading scenarios
            material.rho = 1; % density
            material.mu = 1; % viscosity

            %% Construct FEA Solver
            solver = feaClass(brep,numElements,material, ...
                interpolation,numScenarios); % call superclass

            % inlet
            Uin = 1;
            solver = solver.fixUOfEdge(7,Uin);
            solver = solver.fixVOfEdge(7,0);
            % outlet
            solver = solver.fixPOfEdge(2,0);
            solver = solver.fixUOfEdge(2,0);

            % no-slip top bottom
            solver = solver.fixUOfEdge([1,3,5],0);
            solver = solver.fixVOfEdge([1,3,5],0);
            % no-slip left right
            solver = solver.fixUOfEdge([4,6,8],0);
            solver = solver.fixVOfEdge([4,6,8],0);

            solver = solver.preProcess();

            %% Objective and Constraints
            objective = densityEnergyDissipation(solver);

            volumeFraction = 0.3;
            constraints = {volume(solver, volumeFraction)};

            % manufacturing constraints
            rmin = 1.5;
            mfgConstraints = {
                minimumFeatureSize_dist(solver, rmin)
                physicalDensity(solver)
                };

            %% Construct Optimizer
            testMode = true; % test mode: no export, no plotting
            topopt = topoptClass(solver, ...
                objective,constraints,mfgConstraints, ...
                update, ...
                maxNumIters,exportGIF,testMode);

            %% Optimize
            topopt = topopt.optimize();

            %% Verify convergence
            testCase.verifyTrue(ismember(topopt.m_flag, [1 2]), ...
                'Density TopOpt for thermo-elasticity did not converge to a solution.');
        end

    end
end
