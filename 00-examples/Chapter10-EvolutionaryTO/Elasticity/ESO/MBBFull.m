clc;clear;  close all;format compact; format long
warning('off','all')

%% Solvers
feaClass = @fea2d_elasticity;
topoptClass = @eso2d_elasticity;

%% General Parameters
vectorize = true;
exportImages = false;
exportGIF = false;
exportSTL = false;

%% File Path
p = mfilename("fullpath");
[path,example_name,~] = fileparts(p);

disp("==================================");
disp(['Running ',example_name])

%% Optimizer Parameters
maxNumIters = 300;

%% Problem Definition
brep = 'MBBFull.brep'; % geometry
numElements = 6400; % mesh
material.E = 100e9; material.nu = 0.3; material.rho = 1000; % material
numScenarios = 1;

%% Construct FEA Solver
solver = feaClass(brep,numElements,material,vectorize,numScenarios); % call superclass

solver = solver.fixYOfEdge([1,3]);
solver = solver.fixXOfEdge(1);
solver = solver.applyYForceOnEdge(6,-2e5);

solver = solver.preProcess(); % FEA pre-processing

%% Objective and Constraints
objective = topologicalSensitivityComplianceElasticity(solver);
volumeFraction = 0.5;
constraints  = {volume(solver, volumeFraction)};

% manufacturing constraints
mfgConstraints = {
    minimumFeatureSize_gaussian(solver)
    symmetry_tsf(solver,0) % 0: x-dir, 1: y-dir
    };

%% Construct Optimizer
volDecrement = 0.025;
topopt = topoptClass(solver, ...
    objective,constraints,mfgConstraints, ...
    volDecrement,exportGIF);
%% Make Directory
if exportImages
    folder = [path '/../result/ESO/example' '-' example_name '/']; %#ok
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
ex_title = strjoin({'Evolutionary Topology Optimization for Elasticity ','Example',example_name},' ');
combineFigures(ex_title);
if exportImages
    saveAll(folder);%#ok
end
if exportImages
    diary off
end

cd(path)