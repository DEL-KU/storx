clc;clear;  close all;format compact; format long
warning('off','all')

%% Solvers
feaClass = @fea2d_elasticity;
topoptClass = @density2d_elasticity;

%% General Parameters
vectorize = true;
exportImages = true;
exportGIF = true;

%% File Path
p = mfilename("fullpath");
[path,example_name,~] = fileparts(p);

disp("==================================");
disp(['Running ',example_name])

%% Optimizer Parameters
interpolation = 'simp';
update = 'OC';
maxNumIters = 300;
penaltyStruct = struct('min',3,'max',3,'inc',0.0);

%% Problem Definition
brep = 'GripperComplex.brep'; % geometry
numElements = 10000; % mesh
material.E = 2e9; material.nu = 0.35; material.rho = 1300; % material
force = 10; % N
numScenarios = 1;
%% Construct FEA Solver
solver = feaClass(brep,numElements,material,vectorize,numScenarios, ...
    interpolation,penaltyStruct); % call superclass

solver = solver.fixEdge([5,6,11,12]);
solver = solver.applyXForceOnEdge(18,force);

solver = solver.preProcess(); % FEA pre-processing
solver = solver.solve();
solver = solver.postProcess();
%% Make Directory
if exportImages
    % Make directory
    folder = [path '/../result/example' '-' example_name '/']; %#ok
    mkdir(folder)
    cd(folder)
    delete 'log.txt'
    diary 'log.txt'
end

%% Output
solver.printElascticityResults();
solver.plotGeometryWithLabels();
solver.plotMesh();
solver.plotBoundaryCondition();
solver.plotDeformation();
solver.plotVonMisesStress();
solver.plotPrincipalStress();

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
cd(path)