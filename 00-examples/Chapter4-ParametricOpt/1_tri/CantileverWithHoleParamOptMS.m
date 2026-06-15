clc;clear;  close all;format compact; format long
warning('off','all')

%% General Parameters
exportImages = false;
exportGIF = false;

%% File Path
p = mfilename("fullpath"); 
[path,example_name,~] = fileparts(p);

%% Export
if exportImages
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
params0.value = [0.2 0.15 1.2 0.1];
params0.lb = [0.05 0.05 1 0.05];
params0.ub = [0.4 0.4 1.5 0.2];

objective = 'compliance'; % objective
constraints.area = 1.8; % constraint value
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
    saveAll(folder);%#ok
 end

%% Plot Combined Figures
% ex_title = strjoin({'Parametric Shape Opt. ','Example',example_name},' ');
% combineFigures(ex_title);
% if exportImages
%     saveAll(folder);%#ok
%  end
% cd(path)

if exportImages
    diary off
end
cd(path)

%% Create Problem
function fem = createProblem(brep)
numElements = 1000; % mesh
material.E = 100e9; material.nu = 0.3; material.rho = 1; % material
fem = triFEA2d_elasticity(brep,numElements,material);
fem = fem.fixEdge([2,15]);
fem = fem.applyYForceOnEdge(11,-1e5);
fem = fem.preProcess();
end

%% Create Geometry from Parameters
function geom = createGeom(params)
% params
L = 2.0;%length
H = 1;%height
h = 0.2; % length of force edge

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
