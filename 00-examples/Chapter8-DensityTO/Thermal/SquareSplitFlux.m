clc;clear;  close all;format compact; format long
warning('off','all')

%% Solvers
feaClass = @fea2d_thermal;
topoptClass = @density2d_thermal;

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
topopt = topoptClass(solver, ...
    objective,constraints,mfgConstraints, ...
    update, ...
    maxNumIters,exportGIF);

%% Make Directory
if exportImages
    folder = [path '/../result/example' '-' example_name '/']; %#ok
    name = [update '-' 'numElem' num2str(numElements) '-' 'vf' num2str(volumeFraction)];
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
topopt.m_solver.plotTemperature();

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
ex_title = strjoin({'Desnity TO for Thermal ','Example',example_name},' ');
combineFigures(ex_title);
if exportImages
    saveAll(folder);%#ok
end
cd(path)