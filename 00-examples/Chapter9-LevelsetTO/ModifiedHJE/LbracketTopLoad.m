function topopt = LbracketTopLoad(options)
arguments
    options.vectorize (1,1) logical = true
    options.uniformGrid (1,1) logical = true
    options.exportImages (1,1) logical = false
    options.exportGIF (1,1) logical = false
    options.exportSTL (1,1) logical = false
    options.interpolation (1,:) char = 'none'
    options.maxNumIters (1,1) double {mustBeInteger,mustBePositive} = 1000
    options.penaltyStruct (1,1) struct = struct('min',1,'max',1,'inc',0)
    options.numElements (1,1) double {mustBeInteger,mustBePositive} = 3200
    options.material (1,1) struct = struct('E',100e9,'nu',0.3,'rho',1000)
    options.numScenarios (1,1) double {mustBeInteger,mustBePositive} = 1
    options.load (1,1) double = -1e5
    options.volumeFraction (1,1) double {mustBePositive} = 0.5
    options.topWeight (1,1) double {mustBePositive} = 10
    options.stlThickness (1,1) double {mustBePositive} = 0.2
end

close all; format compact; format long
%% General Parameters
vectorize = options.vectorize;
uniformGrid = options.uniformGrid; % needed for the Hamilton-Jacobi solver
exportImages = options.exportImages;
exportGIF = options.exportGIF;
exportSTL = options.exportSTL;
%% Solvers
feaClass = @fea2d_elasticity;
topoptClass = @modifiedHJ2d_elasticity;

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
brep = 'LBracketNoFillet.brep'; % geometry
numElements = options.numElements; % mesh
material = options.material;
numScenarios = options.numScenarios;

%% Construct FEA Solver
solver = feaClass(brep,numElements,material,vectorize,numScenarios, ...
    interpolation,penaltyStruct,uniformGrid); % call superclass

solver = solver.fixEdge(6);
solver = solver.applyYForceOnEdge(3,options.load);

solver = solver.preProcess(); % FEA pre-processing

%% Objective and Constraints
objective = modifiedHJComplianceElasticity(solver);

volumeFraction = options.volumeFraction;
constraints  = {volume(solver, volumeFraction)};

% manufacturing constraints
mfgConstraints = {minimumFeatureSize_conv(solver)
    retain_levelset(solver,3) }; 

%% Construct Optimizer
topWeight = options.topWeight;
topopt = topoptClass(solver, ...
    objective,constraints,mfgConstraints, ...
    topWeight, ...
    maxNumIters,exportGIF);

%% Make Directory
if exportImages || exportSTL || exportGIF
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
topopt = topopt.optimize();

%% Plotting
topopt.m_solver.plotBoundaryCondition();
topopt.m_solver.plotDeformation();
topopt.m_solver.plotVonMisesStress();
topopt.m_solver.plotPrincipalStress();

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
ex_title = strjoin({'Level-Set Topology Optimization for Elasticity ','Example',example_name},' ');
combineFigures(ex_title);
if exportImages 
    saveAll(folder);%#ok
 end
if exportImages || exportSTL || exportGIF
    diary off
end

cd(path)
end
