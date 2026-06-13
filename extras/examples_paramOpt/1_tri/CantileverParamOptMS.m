clc;clear;  close all;format compact; format long
warning('off','all')

%% General Parameters
exportImages = false;
exportGIF = false;

%% Export
if exportImages
    %% File Path
    p = mfilename("fullpath"); %#ok
    [path,example_name,~] = fileparts(p);
    % Make directory
    folder = [path '/result/example' '-' example_name '/']; %#ok
    mkdir(folder)
    cd(folder)
    delete 'log.txt'
    diary 'log.txt'
end

disp("==================================");
disp(['Running ',example_name])

%% Problem Definition
params0.value = [0.3 0.3];
params0.lb = [0.1 0.1];
params0.ub = [1.5 1.5];

objective = 'compliance'; % objective
constraints.area = 1.4; % constraint value
constraints.type = 'ineq'; % constraint type: 'eq' or 'ineq'
%% Construct Optimizer
brepHandle = @createGeom;
solverHandle = @createProblem;
terminationTolerance = 1e-6;
finiteDifferenceStepSize = 1e-6;

parOpt = parameterOpt2d_MS(brepHandle,solverHandle,params0, ...
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
    saveAll(folder); %#ok
end

%% Plot Combined Figures
ex_title = strjoin({'Parametric Shape Opt. ','Example',example_name},' ');
combineFigures(ex_title);
if exportImages
    saveAll(folder); %#ok
end
cd(path)

%% Create Problem
function fem = createProblem(brep)
numElements = 1000; % mesh
material.E = 100e9; material.nu = 0.3; material.rho = 1; % material
fem = triFEA2d_elasticity(brep,numElements,material);
fem = fem.fixEdge([2,8]);
fem = fem.applyYForceOnEdge(5,-1e5);
end

%% Create Geometry from Parameters
function geom = createGeom(params)
% params
L = 2.0;%length
H = 1;%height
h = 0.2; % length of force edge
hl = 0.2;
a = params(1); % corner cutouts
b = params(2); % left edge cutout
geom.vertices = [b 0;
    0 -hl;
    0 -H/2;
    L-a -H/2;
    L -h/2;
    L h/2;
    L-a H/2;
    0 H/2;
    0 hl ]';
geom.segments = [1 1 2 0; 1 2 3 0;1 3 4 0;1 4 5 0;1 5 6 0;1 6 7 0;1 7 8 0;1 8 9 0;1 9 1 0]';
end
