function shapeopt = SquareSplit(options)
arguments
    options.vectorize (1,1) logical = true
    options.uniformGrid (1,1) logical = true
    options.exportImages (1,1) logical = false
    options.exportGIF (1,1) logical = false
    options.exportSTL (1,1) logical = false
    options.interpolation (1,:) char = 'none'
    options.maxNumIters (1,1) double {mustBeInteger,mustBePositive} = 3000
    options.penaltyStruct (1,1) struct = struct('min',1,'max',1,'inc',0)
    options.numElements (1,1) double {mustBeInteger,mustBePositive} = 5000
    options.material (1,1) struct = struct('k',1)
    options.numScenarios (1,1) double {mustBeInteger,mustBePositive} = 1
    options.flux (1,1) double = 10
    options.volumeFraction (1,1) double {mustBePositive} = 0.5
    options.nHolesX (1,1) double {mustBeInteger,mustBeNonnegative} = 3
    options.nHolesY (1,1) double {mustBeInteger,mustBeNonnegative} = 3
    options.initialHoleRadius (1,1) double {mustBePositive} = 0.5
    options.stlThickness (1,1) double {mustBePositive} = 0.2
end

configureGraphics();

close all;format compact; format long
warning('off','all')

%% Solvers
feaClass = @fea2d_thermal;
shapeoptClass = @standardHJ2d_thermal;

%% General Parameters
vectorize = options.vectorize;
uniformGrid = options.uniformGrid; % needed for the Hamilton-Jacobi solver
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
maxNumIters = options.maxNumIters;
penaltyStruct = options.penaltyStruct;

%% Problem Definition
brep = 'SquareSplit.brep'; % geometry
numElements = options.numElements; % mesh
material = options.material;
numScenarios = options.numScenarios;

%% Construct FEA Solver
solver = feaClass(brep,numElements,material,vectorize,numScenarios, ...
    interpolation,penaltyStruct,uniformGrid); % call superclass

solver = solver.fixEdge(2,0);
solver = solver.applyFlux(5,options.flux);

solver = solver.preProcess();

%% Objective and Constraints
objective = standardHJComplianceThermal(solver);

volumeFraction = options.volumeFraction;
constraints  = {volume(solver, volumeFraction)};

% manufacturing constraints
mfgConstraints = {minimumFeatureSize_conv(solver)
    retain_levelset(solver,5)};

%% Construct Optimizer
nHolesX = options.nHolesX; nHolesY = options.nHolesY; r0 = options.initialHoleRadius;

shapeopt = shapeoptClass(solver, ...
    objective,constraints,mfgConstraints, ...
    nHolesX,nHolesY,r0, ...
    maxNumIters,exportGIF);

%% Make Directory
if exportImages || exportGIF || exportSTL
    folder = [path '/result/example' '-' example_name '/']; %#ok
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
    thickness = options.stlThickness;
    shapeopt.exportSTL(example_name, thickness);
end

%% Plot Combined Figures
ex_title = strjoin({'Level-Set Topology Optimization for Thermal ','Example',example_name},' ');
combineFigures(ex_title);
if exportImages
    saveAll(folder);%#ok
end
if exportImages || exportGIF || exportSTL
    diary off
end

cd(path)
end
