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
maxNumIters = 100;
penaltyStruct = struct('min',3,'max',3,'inc',0);

%% Problem Definition
brep = 'windTunnel.brep'; % geometry
numElements = 40000;             % mesh
numScenarios = 1;               % # loading scenarios

% Specify Reynolds number
Re_in = 10.0;          % desired inlet Reynolds number
Uref  = 1.0;         % reference/inlet velocity

% volume fraction
volumeFraction = 0.85;
activeArea = 1.5*0.5;
Lc = sqrt(volumeFraction*activeArea);

% Non-dimensional material parameters consistent with Re_in
material.rho = 1.0;
material.mu  = material.rho * Uref * Lc / Re_in;

alpha_min = 0.;
alpha_max = material.mu/(1e-5*Lc);
alpha_0 = alpha_min;
qa = 10;
%% Construct FEA Solver
solver = feaClass(brep,numElements,material, ...
    interpolation,numScenarios);

solver = solver.setAlphaValues(alpha_min, alpha_max,alpha_0,qa);

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

%% Pre-Process
solver = solver.preProcess();

%% evaluate reference drag
beta = 1.1;
dragRef = 2.4207; % from Joe Alexandersen's paper
%% Active Design Domain
center = [1.35,0.5];
w = 1.5;
h = 0.5;
solver = solver.createRectangularDesignDomain(center,w,h);
%% Objective and Constraints
objective = densityLift(solver);

constraints = {
    activeVolume(solver, volumeFraction)
    densityDrag(solver,beta*dragRef)
    };

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

% set qa according to the example
topopt = topopt.set_qa(qa);

%% Make Directory
if exportImages
    folder = [path '/../result/lift/example' '-' example_name '/']; %#ok
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
ex_title = strjoin({'Desnity TO for Fluid ','Example',example_name},' ');
combineFigures(ex_title);
if exportImages 
    saveAll(folder);%#ok
 end
if exportImages
    diary off
end

cd(path)