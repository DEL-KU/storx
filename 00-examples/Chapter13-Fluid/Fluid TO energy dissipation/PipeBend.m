clc;clear;  close all;format compact; format long
warning('off','all')

%% Solvers
feaClass = @fea2d_fluid;
topoptClass = @density2d_fluid;

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
update = 'MMA';
maxNumIters = 60;
penaltyStruct = struct('min',3,'max',3,'inc',0);

%% Problem Definition
brep = 'PipeBend.brep'; % geometry
numElements = 2000; % mesh
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
topopt = topoptClass(solver, ...
    objective,constraints,mfgConstraints, ...
    update, ...
    maxNumIters,exportGIF);

%% Make Directory
if exportImages
    folder = [path '/../result/EnergyDissipation/example' '-' example_name '/']; %#ok
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
topopt.m_solver.plotVelocity();
topopt.m_solver.plotPressure();

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
ex_title = strjoin({example_name,'Combined '},' ');
combineFigures(ex_title);
if exportImages
    saveAll(folder);%#ok
end
cd(path)