function topopt = RugbyBall(options)
arguments
    options.exportImages (1,1) logical = false
    options.exportGIF (1,1) logical = false
    options.exportSTL (1,1) logical = false
    options.interpolation (1,:) char = 'simp'
    options.update (1,:) char = 'MMA'
    options.maxNumIters (1,1) double {mustBeInteger,mustBePositive} = 20
    options.penaltyStruct (1,1) struct = struct('min',3,'max',3,'inc',0)
    options.brep (1,:) char = 'Square.brep'
    options.numElements (1,1) double {mustBeInteger,mustBePositive} = 8000
    options.numScenarios (1,1) double {mustBeInteger,mustBePositive} = 1
    options.reynoldsNumber (1,1) double {mustBePositive} = 1000
    options.referenceLength (1,1) double {mustBePositive} = 1
    options.referenceVelocity (1,1) double {mustBePositive} = 1
    options.materialDensity (1,1) double {mustBePositive} = 1
    options.volumeFraction (1,1) double {mustBeGreaterThan(options.volumeFraction,0),mustBeLessThanOrEqual(options.volumeFraction,1)} = 0.85
    options.rmin (1,1) double {mustBePositive} = 1.5
    options.initialCenter (1,2) double = [0.5,0.5]
    options.initialWidth (1,1) double {mustBePositive} = 0.2
    options.initialHeight (1,1) double {mustBePositive} = 0.2
    options.stlThickness (1,1) double {mustBePositive} = 0.2
end
configureGraphics();

close all;format compact; format long
warning('off','all')

%% Solvers
feaClass = @fea2d_fluid;
topoptClass = @density2d_fluid;

%% General Parameters
exportImages = options.exportImages;
exportGIF = options.exportGIF;
exportSTL = options.exportSTL;

%% File Path
p = mfilename("fullpath");
[path,example_name,~] = fileparts(p);

disp("==================================");
disp(['Running ',example_name])

%% Optimizer Parameters
interpolation = options.interpolation;
update = options.update;
maxNumIters = options.maxNumIters;
penaltyStruct = options.penaltyStruct;

%% Problem Definition
brep = options.brep; % geometry
numElements = options.numElements;             % mesh
numScenarios = options.numScenarios;               % # loading scenarios

% Specify Reynolds number
Re_in = options.reynoldsNumber;          % desired inlet Reynolds number
Lref  = options.referenceLength;         % characteristic length (e.g. inlet width), non-dimensional
Uref  = options.referenceVelocity;         % reference/inlet velocity

% Non-dimensional material parameters consistent with Re_in
material.rho = options.materialDensity;
material.mu  = material.rho * Uref * Lref / Re_in;

%% Construct FEA Solver
solver = feaClass(brep,numElements,material, ...
    interpolation,numScenarios,penaltyStruct);

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

volumeFraction = options.volumeFraction;
constraints  = {volume(solver, volumeFraction)};

% manufacturing constraints
rmin = options.rmin;
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
center = options.initialCenter;
w = options.initialWidth;
h = options.initialHeight;
topopt = topopt.setPseudoDensityInRectangle(center,w,h,0,1);

%% Make Directory
if exportImages || exportGIF || exportSTL
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
    thickness = options.stlThickness;
    topopt.exportSTL(example_name, thickness);
end

%% Plot Combined Figures
ex_title = strjoin({example_name,'Combined '},' ');
combineFigures(ex_title);
if exportImages
    saveAll(folder);%#ok
end
if exportImages || exportGIF || exportSTL
    diary off
end

cd(path)
end
