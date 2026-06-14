clear; close all; format compact; format long
%% General Parameters
vectorize = true;
exportImages = false;
exportGIF = false;
exportSTL = false;

%% Solvers
feaClass = @fea2d_elasticity;
topoptClass = @pareto2d_elasticity;

%% File path
p = mfilename("fullpath");
[path,example_name,~] = fileparts(p);

disp("==================================");
disp(['Running ',example_name])

%% Optimizer Parameters
maxNumIters = 300;

%% Problem definition
brep = 'LBracketNoFilletMidLoad.brep'; % geometry
numElements = 3200; % mesh
material.E = 100e9; material.nu = 0.3; material.rho = 1000; % material
numScenarios = 1;
%% Construct FEA Solver
solver = feaClass(brep,numElements,material,vectorize,numScenarios); % call superclass

solver = solver.fixEdge([7,8,9]);
solver = solver.applyYForceOnEdge(3,-1e5);

solver = solver.preProcess(); % FEA pre-processing
%% Objective and Constraints
objective = topologicalSensitivityComplianceElasticity(solver);

volumeFraction = 0.5;
constraints  = {volume(solver, volumeFraction)};

% manufacturing constraints
mfgConstraints = {
    minimumFeatureSize_gaussian(solver)
    };

%% Construct Optimizer
volDecrement = 0.025;
paretoAggressiveness = 0.65;
topopt = topoptClass(solver, ...
    objective,constraints,mfgConstraints, ...
    volDecrement,paretoAggressiveness,exportGIF);

%% Make Directory
if exportImages
    folder = [path '/../result/example' '-' example_name '/']; %#ok
    name = ['numElem' num2str(numElements) '-' 'vf' num2str(volumeFraction)];
    folder = [folder name '/'];
    mkdir(folder)
    cd(folder)
    diary off
    logFile = fullfile(folder, 'log.txt');
    if exist(logFile, 'file')
        delete(logFile)
    end
    diary(logFile)
end

%% Optimize
topopt = topopt.optimize();

%% Plotting
topopt.m_solver.plotBoundaryCondition();
topopt.plotIsoSurface('LS');
topopt.m_solver.plotDeformation();
topopt.m_solver.plotVonMisesStress();
topopt.m_solver.plotPrincipalStress();
topopt.plotConvergence();

%% Save Individual Figures
if exportImages
    saveAll(folder);%#ok
end

%% Export STL
if exportSTL
    thickness = 0.2;
    topopt.exportSTL(example_name, thickness);
end

%% Plot Combined Figures
ex_title = strjoin({'Pareto-tracing Topology Optimization for Elasticity ','Example',example_name},' ');
combineFigures(ex_title);
if exportImages
    saveAll(folder);%#ok
end
if exportImages
    diary off
end

cd(path)