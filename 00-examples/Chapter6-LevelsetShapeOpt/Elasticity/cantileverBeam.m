clc;clear;  close all;format compact; format long
warning('off','all')

%% Solvers
feaClass = @fea2d_elasticity;
shapeoptClass = @standardHJ2d_elasticity;

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
maxNumIters = 500;
penaltyStruct = struct('min',1,'max',1,'inc',0);

%% Problem Definition
brep = 'CantileverBeam.brep'; % geometry
numElements = 3200; % mesh
material.E = 100e9; material.nu = 0.3; material.rho = 1000; % material
numScenarios = 1;

%% Construct FEA Solver
solver = feaClass(brep,numElements,material,vectorize,numScenarios, ...
    interpolation,penaltyStruct,uniformGrid); % call superclass

solver = solver.fixEdge(5);
solver = solver.applyYForceOnEdge(2,-1e5);

solver = solver.preProcess(); % FEA pre-processing

%% Objective and Constraints
objective = standardHJComplianceElasticity(solver);

volumeFraction = 0.5;
constraints  = {volume(solver, volumeFraction)};

% manufacturing constraints
mfgConstraints = {minimumFeatureSize_conv(solver)}; 

%% Construct Optimizer
nHolesX = 0; nHolesY = 0; r0 = 0;

shapeopt = shapeoptClass(solver, ...
    objective,constraints,mfgConstraints, ...
    nHolesX,nHolesY,r0, ...
    maxNumIters,exportGIF);

%% Make Directory
if exportImages || exportSTL || exportGIF
    folder = [path '/../result/example' '-' example_name '/']; %#ok
    name = ['numElem' num2str(numElements) '-' 'vf' num2str(volumeFraction)];
    folder = [folder name '/'];
    mkdir(folder)
    cd(folder)
    delete 'log.txt'
    diary 'log.txt'
end

%% Optimize
shapeopt = shapeopt.optimize();

%% Plotting
shapeopt.m_solver.plotBoundaryCondition();
shapeopt.m_solver.plotDeformation();
shapeopt.m_solver.plotVonMisesStress();
shapeopt.m_solver.plotPrincipalStress();

%% Save Individual Figures
if exportImages 
    saveAll(folder);%#ok
 end

%% Export STL
if exportSTL
    thickness = 0.2;
    shapeopt.exportSTL(example_name, thickness);
end

%% Plot Combined Figures
ex_title = strjoin({'Level-Set Shape Optimization for Elasticity ','Example',example_name},' ');
combineFigures(ex_title);
if exportImages 
    saveAll(folder);%#ok
 end
cd(path)