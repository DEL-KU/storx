clc;clear;  close all;format compact; format long
warning('off','all')

%% Solvers
feaClass = @fea2d_thermal;
shapeoptClass = @standardHJ2d_thermal;

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
maxNumIters = 3000;
penaltyStruct = struct('min',1,'max',1,'inc',0);

%% Problem Definition
brep = 'SquareSplit.brep'; % geometry
numElements = 5000; % mesh
material.k = 1; % material
numScenarios = 1;

%% Construct FEA Solver
solver = feaClass(brep,numElements,material,vectorize,numScenarios, ...
    interpolation,penaltyStruct,uniformGrid); % call superclass

solver = solver.fixEdge(2,0);
solver = solver.applyFlux(5,10);

solver = solver.preProcess();

%% Objective and Constraints
objective = standardHJComplianceThermal(solver);

volumeFraction = 0.5;
constraints = {volume(solver, volumeFraction)};

% manufacturing constraints
mfgConstraints = {minimumFeatureSize_conv(solver)
    retain_levelset(solver,5)
    symmetry_levelset(solver,0) % 0: x-dir, 1: y-dir
    };

%% Construct Optimizer
nHolesX = 0; nHolesY = 0; r0 = 0;

shapeopt = shapeoptClass(solver, ...
    objective,constraints,mfgConstraints, ...
    nHolesX,nHolesY,r0, ...
    maxNumIters,exportGIF);

%% Make Directory
if exportImages
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
shapeopt.m_solver.plotTemperature();

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
ex_title = strjoin({'Level-Set Shape Optimization for Thermal ','Example',example_name},' ');
combineFigures(ex_title);
if exportImages 
    saveAll(folder);%#ok
 end
cd(path)