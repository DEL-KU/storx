function parOpt = CantileverWithHoleParamOptGS(options)
arguments
    options.exportImages (1,1) logical = false
    options.exportGIF (1,1) logical = false
    options.params0 (1,1) struct = struct('value',[0.1 0.15 1.2 0.1], ...
        'lb',[0.05 0.05 1 0.05],'ub',[0.4 0.4 1.5 0.2])
    options.objective (1,:) char = 'compliance'
    options.constraints (1,1) struct = struct('area',1.8,'type','ineq')
    options.terminationTolerance (1,1) double {mustBePositive} = 1e-6
    options.finiteDifferenceStepSize (1,1) double {mustBePositive} = 1e-6
    options.vectorize (1,1) logical = true
    options.numElements (1,1) double {mustBeInteger,mustBePositive} = 1000
    options.material (1,1) struct = struct('E',100e9,'nu',0.3,'rho',1)
    options.numScenarios (1,1) double {mustBeInteger,mustBePositive} = 1
    options.fixedEdges (1,:) double {mustBeInteger,mustBePositive} = [2 15]
    options.forceEdge (1,1) double {mustBeInteger,mustBePositive} = 11
    options.force (1,1) double = -1e5
    options.geometryLength (1,1) double {mustBePositive} = 2.0
    options.geometryHeight (1,1) double {mustBePositive} = 1
    options.loadEdgeLength (1,1) double {mustBePositive} = 0.2
end

configureGraphics();
close all;format compact; format long
warning('off','all')

%% General Parameters
exportImages = options.exportImages;
exportGIF = options.exportGIF;

%% File Path
p = mfilename("fullpath"); 
[path,example_name,~] = fileparts(p);

%% Export
if exportImages || exportGIF
    % Make directory
    folder = [path '/result/example' '-' example_name '/']; %#ok
    mkdir(folder)
    cd(folder)
    diary off
    logFile = fullfile(folder, 'log.txt');
    if exist(logFile, 'file')
        delete(logFile)
    end
    diary(logFile)
end

disp("==================================");
disp(['Running ',example_name])

%% Problem Definition
params0 = options.params0;

objective = options.objective; % objective
constraints = options.constraints;
%% Construct Optimizer
brepHandle = @(params) createGeom(params,options.geometryLength, ...
    options.geometryHeight,options.loadEdgeLength);
solverHandle = @(brep) createProblem(brep,options.vectorize, ...
    options.numElements,options.material,options.numScenarios, ...
    options.fixedEdges,options.forceEdge,options.force);
terminationTolerance = options.terminationTolerance;
finiteDifferenceStepSize = options.finiteDifferenceStepSize;

parOpt = parameterOpt2d_GS(brepHandle,solverHandle,params0, ...
    objective,constraints, ...
    terminationTolerance,finiteDifferenceStepSize,exportGIF);

%% Optimize
parOpt = parOpt.optimize();
%% Output
parOpt.m_solverInitial.plotGeometry(1,0, 'Initial Geometry');
parOpt.m_solverFinal.plotDeformation();
parOpt.m_solverFinal.plotVonMisesStress();
%% Save
if exportImages 
    saveAll(folder);%#ok
 end

%% Plot Combined Figures
ex_title = strjoin({'Parametric Shape Opt. ','Example',example_name},' ');
combineFigures(ex_title);
if exportImages 
    saveAll(folder);%#ok
 end
if exportImages || exportGIF
    diary off
end

cd(path)
end

%% Create Problem
function fem = createProblem(brep,vectorize,numElements,material,numScenarios, ...
        fixedEdges,forceEdge,force)
fem = fea2d_elasticity(brep,numElements,material,vectorize,numScenarios);
fem = fem.fixEdge(fixedEdges);
fem = fem.applyYForceOnEdge(forceEdge,force);
fem = fem.preProcess();
end

%% Create Geometry from Parameters
function geom = createGeom(params,L,H,h)
a = params(1); % corner cutouts
b = params(2); % left edge cutout
c = params(3);
r = params(4);
geom.vertices = [b 0;
    0 -b; 
    0 -H/2;
    c -H/2;
    c -r;
    c r;
    c 0;
    L-a -H/2;
    L -H/2+a; 
    L -h/2; 
    L h/2; 
    L H/2-a; 
    L-a H/2; 
    0 H/2;
    0 b ]';

geom.segments = [1 1 2 0; 
    1 2 3 0;
    1 3 4 0;
    -1 4 5 0
    2 5 6 7;
    2 6 5 7;
    -1 5 4 0;
    1 4 8 0;
    1 8 9 0;
    1 9 10 0;
    1 10 11 0;
    1 11 12 0;
    1 12 13 0
    1 13 14 0
    1 14 15 0
    1 15 1 0]';
end
