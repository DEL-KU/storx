function topopt = windTunnel_drag(options)
arguments
    options.exportImages (1,1) logical = false
    options.exportGIF (1,1) logical = false
    options.interpolation (1,:) char = 'simp'
    options.update (1,:) char = 'MMA'
    options.maxNumIters (1,1) double {mustBeInteger,mustBePositive} = 100
    options.penaltyStruct (1,1) struct = struct('min',3,'max',3,'inc',0)
    options.brep (1,:) char = 'windTunnel.brep'
    options.numElements (1,1) double {mustBeInteger,mustBePositive} = 40000
    options.numScenarios (1,1) double {mustBeInteger,mustBePositive} = 1
    options.reynoldsNumber (1,1) double {mustBePositive} = 10
    options.referenceVelocity (1,1) double {mustBePositive} = 1
    options.volumeFraction (1,1) double {mustBeGreaterThan(options.volumeFraction,0),mustBeLessThanOrEqual(options.volumeFraction,1)} = 0.85
    options.materialDensity (1,1) double {mustBePositive} = 1
    options.alphaMin (1,1) double {mustBeGreaterThanOrEqual(options.alphaMin,0)} = 0
    options.alphaScale (1,1) double {mustBePositive} = 1e-5
    options.qa (1,1) double {mustBePositive} = 10
    options.designCenter (1,2) double = [1.35,0.5]
    options.designWidth (1,1) double {mustBePositive} = 1.5
    options.designHeight (1,1) double {mustBePositive} = 0.5
    options.rmin (1,1) double {mustBePositive} = 1.5
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
Uref  = options.referenceVelocity;         % reference/inlet velocity

% volume fraction
volumeFraction = options.volumeFraction;
activeArea = options.designWidth*options.designHeight;
Lc = sqrt(volumeFraction*activeArea);

% Non-dimensional material parameters consistent with Re_in
material.rho = options.materialDensity;
material.mu  = material.rho * Uref * Lc / Re_in;

alpha_min = options.alphaMin;
alpha_max = material.mu/(options.alphaScale*Lc);
alpha_0 = alpha_min;
qa = options.qa;
%% Construct FEA Solver
solver = feaClass(brep,numElements,material, ...
    interpolation,numScenarios,penaltyStruct);

solver = solver.setAlphaValues(alpha_min, alpha_max,alpha_0,qa);

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

%% Pre-Process
solver = solver.preProcess();

%% Active Design Domain
center = options.designCenter;
w = options.designWidth;
h = options.designHeight;
solver = solver.createRectangularDesignDomain(center,w,h);
%% Objective and Constraints
objective = densityDrag(solver);

constraints = {
    activeVolume(solver, volumeFraction)
    };

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

% set qa according to the example
topopt = topopt.set_qa(qa);
%% Make Directory
if exportImages || exportGIF
    folder = [path '/../result/drag/example' '-' example_name '/']; %#ok
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

%% Plot Combined Figures
ex_title = strjoin({'Desnity TO for Fluid ','Example',example_name},' ');
combineFigures(ex_title);
if exportImages 
    saveAll(folder);%#ok
 end
if exportImages || exportGIF
    diary off
end

cd(path)
end
