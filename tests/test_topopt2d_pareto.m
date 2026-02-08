classdef test_topopt2d_pareto < matlab.unittest.TestCase
    %% Test class for Pareto-tracing topopt
    methods (Test)

        function testLbrackt_pareto(testCase)
            testMode = true;
            %% General Parameters
            vectorize = true;
            exportGIF = false;

            %% Solvers
            feaClass = @fea2d_elasticity;
            topoptClass = @pareto2d_elasticity;

            %% Problem definition
            brep = 'LBracketNoFillet.brep'; % geometry
            numElements = 3200; % mesh
            material.E = 100e9; material.nu = 0.3; material.rho = 1; % material
            numScenarios = 1;

            %% Construct FEA Solver
            solver = feaClass(brep,numElements,material,vectorize,numScenarios); % call superclass

            solver = solver.fixEdge(6);
            solver = solver.applyYForceOnEdge(3,-1e5);

            solver = solver.preProcess(); % FEA pre-processing

            %% Objective and Constraints
            objective = topologicalSensitivityComplianceElasticity(solver);

            volumeFraction = 0.5;
            constraints = {volume(solver, volumeFraction)};

            % manufacturing constraints
            mfgConstraints = {
                minimumFeatureSize_gaussian(solver)
                };

            %% Construct Optimizer
            volDecrement = 0.025;
            paretoAggressiveness = 0.65;
            topopt = topoptClass(solver, ...
                objective,constraints,mfgConstraints, ...
                volDecrement,paretoAggressiveness,exportGIF,testMode);

            %% Optimize
            topopt = topopt.optimize();

            %% Verify convergence
            testCase.verifyTrue(ismember(topopt.m_flag, [1 2]), ...
                'Pareto-tracing Topology Optimization for elasticity did not converge to a solution.');

        end

    end
end
