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
maxNumIters = 20;
penaltyStruct = struct('min',3,'max',3,'inc',0);

%% Problem Definition
brep = 'Square.brep'; % geometry
numElements = 8000;             % mesh
numScenarios = 1;               % # loading scenarios

% Specify Reynolds number
Re_in = 1000;          % desired inlet Reynolds number
Lref  = 1.0;         % characteristic length (e.g. inlet width), non-dimensional
Uref  = 1.0;         % reference/inlet velocity

% Non-dimensional material parameters consistent with Re_in
material.rho = 1.0;
material.mu  = material.rho * Uref * Lref / Re_in;

%% Construct FEA Solver
solver = feaClass(brep,numElements,material, ...
    interpolation,numScenarios);

% inlet
Uin = Uref;  % controlled indirectly by Re_in via mu
solver = solver.fixUOfEdge([4,5,6],Uin,0); % uniform profile
solver = solver.fixVOfEdge([4,5,6],0,0);
% outlet
solver = solver.fixPOfEdge([4,6],0);
% no-slip top bottom
solver = solver.fixVOfEdge([1,3],0);
% no-slip left right
solver = solver.fixUOfEdge(2,Uin,0); % uniform profile
solver = solver.fixVOfEdge(2,0,0);

solver = solver.preProcess();
%% Objective and Constraints
objective = densityEnergyDissipation(solver);

volumeFraction = 0.85;
constraints  = {volume(solver, volumeFraction)};

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

% set initial solid
center = [0.5,0.5];
w = 0.2;
h = 0.2;
topopt = topopt.setPseudoDensityInRectangle(center,w,h,0,1);

%% Make Directory
if exportImages
    folder = [path '/../result/EnergyDissipation/example' '-' example_name '/']; %#ok
    name = [update '-' 'numElem' num2str(numElements) '-' 'vf' num2str(volumeFraction)];
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
if exportImages
    diary off
end

cd(path)