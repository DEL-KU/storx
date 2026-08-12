function topopt = MBBSymmetry(options)
arguments
    options.vectorize (1,1) logical = true
    options.exportImages (1,1) logical = false
    options.combineFigure (1,1) logical = false
    options.exportGIF (1,1) logical = false
    options.exportSTL (1,1) logical = false
    options.interpolation (1,:) char = 'simp'
    options.update (1,:) char = 'OC'
    options.maxNumIters (1,1) double {mustBeInteger,mustBePositive} = 300
    options.penaltyStruct (1,1) struct = struct('min',3,'max',3,'inc',0)
    options.numElements (1,1) double {mustBeInteger,mustBePositive} = 3200
    options.material (1,1) struct = struct('E',100e9,'nu',0.3,'rho',1000)
    options.numScenarios (1,1) double {mustBeInteger,mustBePositive} = 1
    options.load (1,1) double = -1e5
    options.volumeFraction (1,1) double {mustBePositive} = 0.5
    options.filterRadius (1,1) double {mustBePositive} = 1.5
    options.stlThickness (1,1) double {mustBePositive} = 0.1
end

configureGraphics();

close all;format compact; format long
warning('off','all')

%% Solvers
feaClass = @fea2d_elasticity;
topoptClass = @density2d_elasticity;

%% General Parameters
vectorize = options.vectorize;
exportImages = options.exportImages;
combineFigure = options.combineFigure;
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
brep = 'MBBSymmetry.brep'; % geometry
numElements = options.numElements; % mesh
material = options.material;
numScenarios = options.numScenarios;

%% Construct FEA Solver
solver = feaClass(brep,numElements,material,vectorize,numScenarios, ...
    interpolation,penaltyStruct); % call superclass

solver = solver.fixYOfEdge(2);
solver = solver.fixXOfEdge(6);
solver = solver.applyYForceOnEdge(5,options.load);

solver = solver.preProcess(); % FEA pre-processing

%% Objective and Constraints
objective = densityComplianceElasticity(solver);

volumeFraction = options.volumeFraction;
constraints  = {volume(solver, volumeFraction)};

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
if exportImages || exportGIF || exportSTL
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
if combineFigure
    ex_title = strjoin({example_name,'Combined '},' ');
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
