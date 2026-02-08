clc;clear;  close all;format compact; format long
warning('off','all')

%% Solvers
feaClass = @fea2d_elasticity;
topoptClass = @modifiedHJ2d_elasticity;

%% General Parameters
vectorize = true;
uniformGrid = 1; % needed for the Hamilton-Jacobi solver
exportImages = false;
exportGIF = false;
exportSTL = false;
%% File Path
p = mfilename("fullpath");
[path,example_name,~] = fileparts(p);

disp("==================================");
disp(['Running ',example_name])

%% Optimizer Parameters
interpolation = 'none';
maxNumIters = 300;
penaltyStruct = struct('min',1,'max',1,'inc',0);

%% Problem Definition
brep = 'MBBSymmetry.brep'; % geometry
numElements = 3200; % mesh
material.E = 100e9; material.nu = 0.3; material.rho = 1000; % material
numScenarios = 1;

%% Construct FEA Solver
solver = feaClass(brep,numElements,material,vectorize,numScenarios, ...
    interpolation,penaltyStruct,uniformGrid); % call superclass

solver = solver.fixYOfEdge(2);
solver = solver.fixXOfEdge(6);
solver = solver.applyYForceOnEdge(5,-1e5);

solver = solver.preProcess(); % FEA pre-processing
%% Objective and Constraints
objective = modifiedHJComplianceElasticity(solver);

volumeFraction = 0.5;
constraints  = {volume(solver, volumeFraction)};

% manufacturing constraints
mfgConstraints = {
    minimumFeatureSize_conv(solver,1)
    retain_levelset(solver,5)
    };

%% Construct Optimizer
topWeight = 10;
topopt = topoptClass(solver, ...
    objective,constraints,mfgConstraints, ...
    topWeight, ...
    maxNumIters,exportGIF);

%% Make Directory
if exportImages
    folder = [path '/result/example' '-' example_name '/']; %#ok
    name = ['numElem' num2str(numElements) '-' 'vf' num2str(volumeFraction)];
    folder = [folder name '/'];
    mkdir(folder)
    cd(folder)
    delete 'log.txt'
    diary 'log.txt'
end

%% Optimize
topopt = topopt.optimize();

%% Plotting
topopt.m_solver.plotBoundaryCondition();
topopt.m_solver.plotDeformation();
topopt.m_solver.plotVonMisesStress();
topopt.m_solver.plotPrincipalStress();

%% Save Individual Figures
if exportImages
    saveAll(folder);%#ok
end

%% Export STL
if exportSTL
    thickness = 0.2;
    topopt.exportSTL(folder,example_name, thickness);
end

%% Plot Combined Figures
ex_title = strjoin({'Level-Set Topology Optimization for Elasticity ','Example',example_name},' ');
combineFigures(ex_title);
if exportImages
    saveAll(folder);%#ok
end
cd(path)