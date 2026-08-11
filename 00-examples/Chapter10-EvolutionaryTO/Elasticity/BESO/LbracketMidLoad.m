function topopt = LbracketMidLoad(options)
arguments
    options.vectorize (1,1) logical = true
    options.exportImages (1,1) logical = false
    options.exportGIF (1,1) logical = false
    options.exportSTL (1,1) logical = false
    options.brep (1,:) char = 'LBracketNoFilletMidLoad.brep'
    options.numElements (1,1) double {mustBeInteger,mustBePositive} = 3200
    options.material (1,1) struct = struct('E',100e9,'nu',0.3,'rho',1000)
    options.numScenarios (1,1) double {mustBeInteger,mustBePositive} = 1
    options.loadMagnitude (1,1) double = -1e5
    options.volumeFraction (1,1) double {mustBeGreaterThan(options.volumeFraction,0),mustBeLessThanOrEqual(options.volumeFraction,1)} = 0.5
    options.volDecrement (1,1) double {mustBeGreaterThan(options.volDecrement,0),mustBeLessThanOrEqual(options.volDecrement,1)} = 0.025
    options.filterSigma (1,1) double {mustBePositive} = 0.6
    options.stlThickness (1,1) double {mustBePositive} = 0.2
end
configureGraphics();

close all; format compact; format long
%% General Parameters
vectorize = options.vectorize;
exportImages = options.exportImages;
exportGIF = options.exportGIF;
exportSTL = options.exportSTL;

%% Solvers
feaClass = @fea2d_elasticity;
topoptClass = @beso2d_elasticity;

%% File path
p = mfilename("fullpath");
[path,example_name,~] = fileparts(p);

disp("==================================");
disp(['Running ',example_name])

%% Optimizer Parameters

%% Problem definition
brep = options.brep; % geometry
numElements = options.numElements; % mesh
material = options.material; % material
numScenarios = options.numScenarios;

%% Construct FEA Solver
solver = feaClass(brep,numElements,material,vectorize,numScenarios); % call superclass

solver = solver.fixEdge([7,8,9]);
solver = solver.applyYForceOnEdge(3,options.loadMagnitude);

solver = solver.preProcess(); % FEA pre-processing
%% Objective and Constraints
objective = topologicalSensitivityComplianceElasticity(solver);
volumeFraction = options.volumeFraction;
constraints  = {volume(solver, volumeFraction)};

% manufacturing constraints
mfgConstraints = {
    minimumFeatureSize_gaussian(solver,options.filterSigma)
    };

%% Construct Optimizer
volDecrement = options.volDecrement;
topopt = topoptClass(solver, ...
    objective,constraints,mfgConstraints, ...
    volDecrement,exportGIF);

%% Make Directory
if exportImages || exportGIF || exportSTL
    folder = [path '/../result/BESO/example' '-' example_name '/']; %#ok
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
ex_title = strjoin({'Evolutionary Topology Optimization for Elasticity ','Example',example_name},' ');
combineFigures(ex_title);
if exportImages
    saveAll(folder);%#ok
end
if exportImages || exportGIF || exportSTL
    diary off
end

cd(path)
end
