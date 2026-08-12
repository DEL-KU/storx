function shapeopt = LbracketMidLoad(options)
arguments
    options.vectorize (1,1) logical = true
    options.uniformGrid (1,1) double {mustBeInteger,mustBePositive} = 1
    options.exportImages (1,1) logical = false
    options.combineFigure (1,1) logical = false
    options.exportGIF (1,1) logical = false
    options.exportSTL (1,1) logical = false
    options.interpolation (1,:) char = 'none'
    options.maxNumIters (1,1) double {mustBeInteger,mustBePositive} = 3200
    options.penaltyStruct (1,1) struct = struct('min',1,'max',1,'inc',0)
    options.brep (1,:) char = 'LBracketNoFilletMidLoad.brep'
    options.numElements (1,1) double {mustBeInteger,mustBePositive} = 3200
    options.material (1,1) struct = struct('E',100e9,'nu',0.3,'rho',1000)
    options.numScenarios (1,1) double {mustBeInteger,mustBePositive} = 1
    options.volumeFraction (1,1) double {mustBePositive} = 0.5
    options.nHolesX (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    options.nHolesY (1,1) double {mustBeInteger,mustBeNonnegative} = 0
    options.r0 (1,1) double {mustBeNonnegative} = 0
    options.stlThickness (1,1) double {mustBePositive} = 0.2
end

configureGraphics();

close all; format compact; format long
%% General Parameters
vectorize = options.vectorize;
uniformGrid = options.uniformGrid; % needed for the Hamilton-Jacobi solver
exportImages = options.exportImages;
combineFigure = options.combineFigure;
exportGIF = options.exportGIF;
exportSTL = options.exportSTL;

%% Solvers
feaClass = @fea2d_elasticity;
shapeoptClass = @standardHJ2d_elasticity;

%% File path
p = mfilename("fullpath"); 
[path,example_name,~] = fileparts(p);

disp("==================================");
disp(['Running ',example_name])

%% Optimizer Parameters
interpolation = options.interpolation;
maxNumIters = options.maxNumIters;
penaltyStruct = options.penaltyStruct;

%% Problem definition
brep = options.brep; % geometry
numElements = options.numElements; % mesh
material = options.material; % material
numScenarios = options.numScenarios;

%% Construct FEA Solver
solver = feaClass(brep,numElements,material,vectorize,numScenarios, ...
    interpolation,penaltyStruct,uniformGrid); % call superclass

solver = solver.fixEdge([7,8,9]);
solver = solver.applyYForceOnEdge(3,-1e5);

solver = solver.preProcess(); % FEA pre-processing

%% Objective and Constraints
objective = standardHJComplianceElasticity(solver);

volumeFraction = options.volumeFraction;
constraints  = {volume(solver, volumeFraction)};

% manufacturing constraints
mfgConstraints = {minimumFeatureSize_conv(solver)
    retain_levelset(solver,[3])};

%% Construct Optimizer
nHolesX = options.nHolesX; nHolesY = options.nHolesY; r0 = options.r0;

shapeopt = shapeoptClass(solver, ...
    objective,constraints,mfgConstraints, ...
    nHolesX,nHolesY,r0, ...
    maxNumIters,exportGIF);

%% Make Directory
if exportImages || exportGIF || exportSTL
    folder = [path '/../result/example' '-' example_name '/']; %#ok
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
shapeopt.m_solver.plotDeformation();
shapeopt.m_solver.plotVonMisesStress();
shapeopt.m_solver.plotPrincipalStress();

%% Save Individual Figures
if exportImages 
    saveAll(folder);%#ok
 end

 %% Export STL
if exportSTL
    shapeopt.exportSTL(example_name, options.stlThickness);
end

%% Plot Combined Figures
if combineFigure
    ex_title = strjoin({'Level-Set Shape Optimization for Elasticity ','Example',example_name},' ');
    combineFigures(ex_title);
    if exportImages
        saveAll(folder);%#ok
    end
end
if exportImages || exportGIF || exportSTL
    diary off
end

cd(path)
end
