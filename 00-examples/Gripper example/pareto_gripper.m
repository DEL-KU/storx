function topopt = pareto_gripper(options)
arguments
    options.vectorize (1,1) logical = true
    options.exportImages (1,1) logical = false
    options.exportGIF (1,1) logical = false
    options.exportSTL (1,1) logical = false
    options.numElements (1,1) double {mustBeInteger,mustBePositive} = 10000
    options.material (1,1) struct = struct('E',2e9,'nu',0.35,'rho',1300)
    options.force (1,1) double = 10
    options.numScenarios (1,1) double {mustBeInteger,mustBePositive} = 1
    options.volumeFraction (1,1) double {mustBePositive} = 0.65
    options.volumeDecrement (1,1) double {mustBePositive} = 0.025
    options.paretoAggressiveness (1,1) double {mustBePositive} = 0.8
    options.stlThickness (1,1) double {mustBePositive} = 10
    options.stlMinPoints (1,1) double {mustBeInteger,mustBeNonnegative} = 10
end

clc; close all;format compact; format long
warning('off','all')

%% Solvers
feaClass = @fea2d_elasticity;
topoptClass = @pareto2d_elasticity;

%% General Parameters
vectorize = options.vectorize;
exportImages = options.exportImages;
exportGIF = options.exportGIF;
exportSTL = options.exportSTL;

%% File Path
p = mfilename("fullpath");
[path,example_name,~] = fileparts(p);

disp("==================================");
disp(['Running ',example_name])

%% Problem Definition
brep = 'GripperComplex.brep'; % geometry
numElements = options.numElements; % mesh
material = options.material;
force = options.force; % N
numScenarios = options.numScenarios;
%% Construct FEA Solver
solver = feaClass(brep,numElements,material,vectorize,numScenarios); % call superclass

solver = solver.fixEdge([5,6,11,12]);
solver = solver.applyXForceOnEdge(18,force);

solver = solver.preProcess(); % FEA pre-processing

%% Objective and Constraints
objective = topologicalSensitivityComplianceElasticity(solver);

volumeFraction = options.volumeFraction;
constraints  = {volume(solver, volumeFraction)};

% manufacturing constraints
mfgConstraints = {
    minimumFeatureSize_gaussian(solver)
    retain_tsf(solver,[5,6,11,12,18])
    }; 

%% Construct Optimizer
volDecrement = options.volumeDecrement;
paretoAggressiveness = options.paretoAggressiveness;
topopt = topoptClass(solver, ...
    objective,constraints,mfgConstraints, ...
    volDecrement,paretoAggressiveness,exportGIF);

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
topopt.plotIsoSurface('Contour');
topopt.m_solver.plotDeformation();
topopt.m_solver.plotVonMisesStress();
topopt.m_solver.plotPrincipalStress();
topopt.plotConvergence();
%% Export STL
if exportSTL
    thickness = options.stlThickness;
    minPts = options.stlMinPoints;
    topopt.exportSTL(example_name, thickness,minPts);
end

%% Save Individual Figures
if exportImages
    saveAll(folder);%#ok
end

%% Plot Combined Figures
ex_title = strjoin({example_name,'Combined '},' ');
combineFigures(ex_title);
if exportImages
    saveAll(folder);%#ok
end
if exportImages || exportSTL || exportGIF
    diary off
end

cd(path)
end
