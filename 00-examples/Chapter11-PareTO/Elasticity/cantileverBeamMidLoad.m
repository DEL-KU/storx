function topopt = cantileverBeamMidLoad(options)
arguments
    options.vectorize (1,1) logical = true
    options.exportImages (1,1) logical = false
    options.exportGIF (1,1) logical = false
    options.exportSTL (1,1) logical = false
    options.brep (1,:) char = 'CantileverBeamMid.brep'
    options.numElements (1,1) double {mustBeInteger,mustBePositive} = 3200
    options.material (1,1) struct = struct('E',100e9,'nu',0.3,'rho',1000)
    options.numScenarios (1,1) double {mustBeInteger,mustBePositive} = 1
    options.loadMagnitude (1,1) double = -1e5
    options.volumeFraction (1,1) double {mustBeGreaterThan(options.volumeFraction,0),mustBeLessThanOrEqual(options.volumeFraction,1)} = 0.5
    options.volDecrement (1,1) double {mustBeGreaterThan(options.volDecrement,0),mustBeLessThanOrEqual(options.volDecrement,1)} = 0.025
    options.paretoAggressiveness (1,1) double {mustBeGreaterThan(options.paretoAggressiveness,0),mustBeLessThanOrEqual(options.paretoAggressiveness,1)} = 0.65
    options.filterSigma (1,1) double {mustBePositive} = 0.6
    options.stlThickness (1,1) double {mustBePositive} = 0.2
end
configureGraphics();

close all;format compact; format long
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

%% Optimizer Parameters

%% Problem Definition
brep = options.brep; % geometry
numElements = options.numElements; % mesh
material = options.material; % material
numScenarios = options.numScenarios;

%% Construct FEA Solver
solver = feaClass(brep,numElements,material,vectorize,numScenarios); % call superclass

solver = solver.fixEdge(6);
solver = solver.applyYForceOnEdge(3,options.loadMagnitude);

solver = solver.preProcess(); % FEA pre-processing

%% Objective and Constraints
objective = topologicalSensitivityComplianceElasticity(solver);

volumeFraction = options.volumeFraction;
constraints  = {volume(solver, volumeFraction)};

% manufacturing constraints
mfgConstraints = {
    minimumFeatureSize_gaussian(solver,options.filterSigma)
    symmetry_tsf(solver,1) % 0: x-dir, 1: y-dir
    };

%% Construct Optimizer
volDecrement = options.volDecrement;
paretoAggressiveness = options.paretoAggressiveness;
topopt = topoptClass(solver, ...
    objective,constraints,mfgConstraints, ...
    volDecrement,paretoAggressiveness,exportGIF);

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
topopt = topopt.optimize();

%% Plotting
topopt.m_solver.plotBoundaryCondition();
topopt.plotIsoSurface('LS');
topopt.m_solver.plotDeformation();
topopt.m_solver.plotVonMisesStress();
topopt.m_solver.plotPrincipalStress();
topopt.plotConvergence();

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
ex_title = strjoin({'Pareto-tracing Topology Optimization for Elasticity ','Example',example_name},' ');
combineFigures(ex_title);
if exportImages
    saveAll(folder);%#ok
end
if exportImages || exportGIF || exportSTL
    diary off
end

cd(path)
end
