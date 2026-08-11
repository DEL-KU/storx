function topopt = DoublePipe(options)
arguments
    options.exportImages (1,1) logical = false
    options.exportGIF (1,1) logical = false
    options.exportSTL (1,1) logical = false
    options.interpolation (1,:) char = 'simp'
    options.update (1,:) char = 'MMA'
    options.maxNumIters (1,1) double {mustBeInteger,mustBePositive} = 10
    options.penaltyStruct (1,1) struct = struct('min',3,'max',3,'inc',0)
    options.brep (1,:) char = 'DoublePipe.brep'
    options.numElements (1,1) double {mustBeInteger,mustBePositive} = 2000
    options.numScenarios (1,1) double {mustBeInteger,mustBePositive} = 1
    options.material (1,1) struct = struct('rho',1,'mu',1)
    options.inletVelocity (1,1) double = 1
    options.volumeFraction (1,1) double {mustBeGreaterThan(options.volumeFraction,0),mustBeLessThanOrEqual(options.volumeFraction,1)} = 0.3
    options.rmin (1,1) double {mustBePositive} = 1.5
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
numElements = options.numElements; % mesh
numScenarios = options.numScenarios; % # loading scenarios
material = options.material;

%% Construct FEA Solver
solver = feaClass(brep,numElements,material, ...
    interpolation,numScenarios,penaltyStruct); % call superclass

% inlet
Uin = options.inletVelocity;
solver = solver.fixUOfEdge([9,11],Uin);
solver = solver.fixVOfEdge([9,11],0);
% outlet
solver = solver.fixPOfEdge([3,5],0);
solver = solver.fixVOfEdge([3,5],0);

% no-slip top bottom
solver = solver.fixUOfEdge([1,7],0);
solver = solver.fixVOfEdge([1,7],0);
% no-slip left right
solver = solver.fixUOfEdge([2,4,6,8,10,12],0);
solver = solver.fixVOfEdge([2,4,6,8,10,12],0);

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

%% Make Directory
if exportImages || exportGIF || exportSTL
    folder = [path '/../result/EnergyDiss/example' '-' example_name '/']; %#ok
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
