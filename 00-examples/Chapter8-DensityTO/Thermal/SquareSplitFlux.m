function topopt = SquareSplitFlux(options)
arguments
    options.vectorize (1,1) logical = true
    options.exportImages (1,1) logical = false
    options.exportGIF (1,1) logical = false
    options.exportSTL (1,1) logical = false
    options.interpolation (1,:) char = 'simp'
    options.update (1,:) char = 'OC'
    options.maxNumIters (1,1) double {mustBeInteger,mustBePositive} = 500
    options.penaltyStruct (1,1) struct = struct('min',3,'max',3,'inc',0)
    options.numElements (1,1) double {mustBeInteger,mustBePositive} = 10000
    options.material (1,1) struct = struct('k',1)
    options.numScenarios (1,1) double {mustBeInteger,mustBePositive} = 1
    options.flux (1,1) double = 1
    options.volumeFraction (1,1) double {mustBePositive} = 0.5
    options.filterRadius (1,1) double {mustBePositive} = 1.5
    options.stlThickness (1,1) double {mustBePositive} = 0.2
end

configureGraphics();

close all;format compact; format long
warning('off','all')

%% Solvers
feaClass = @fea2d_thermal;
topoptClass = @density2d_thermal;

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
interpolation = options.interpolation;
update = options.update;
maxNumIters = options.maxNumIters;
penaltyStruct = options.penaltyStruct;

%% Problem Definition
brep = 'SquareSplit.brep'; % geometry
numElements = options.numElements; % mesh
numScenarios = options.numScenarios; % # loading scenarios
material = options.material;
%% Construct FEA Solver
solver = feaClass(brep,numElements,material,vectorize,numScenarios, ...
    interpolation,penaltyStruct); % call superclass

solver = solver.fixEdge(2,0);
solver = solver.applyFlux(5,options.flux);

solver = solver.preProcess(); % FEA pre-processing

%% Objective and Constraints
objective = densityComplianceThermal(solver);

volumeFraction = options.volumeFraction;
constraints = {volume(solver, volumeFraction)};

% manufacturing constraints
rmin = options.filterRadius;
mfgConstraints = {
    minimumFeatureSize_dist(solver, rmin)
    };

%% Construct Optimizer
topopt = topoptClass(solver, ...
    objective,constraints,mfgConstraints, ...
    update, ...
    maxNumIters,exportGIF);

%% Make Directory
if exportImages || exportSTL || exportGIF
    folder = [path '/../result/example' '-' example_name '/']; %#ok
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
topopt.m_solver.plotTemperature();

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
ex_title = strjoin({'Desnity TO for Thermal ','Example',example_name},' ');
combineFigures(ex_title);
if exportImages
    saveAll(folder);%#ok
end
if exportImages || exportSTL || exportGIF
    diary off
end

cd(path)
end
