clc;clear;  close all;format compact; format long
warning('off','all')

%% Solvers
feaClass = @fea2d_elasticity;
topoptClass = @density2d_elasticity;

%% General Parameters
vectorize = true;
exportImages = false;
exportGIF = false;
exportSTL = false;
%% File Path
p = mfilename("fullpath");
[path,example_name,~] = fileparts(p);

disp("==================================");
disp(['Running ',example_name])

%% Optimizer Parameters
interpolation = 'simp';
update = 'OC';
maxNumIters = 300;
penaltyStruct = struct('min',3,'max',3,'inc',0);

%% Problem Definition
brep = 'CantileverBeamMid.brep'; % geometry
numElements = 3200; % mesh
material.E = 100e9; material.nu = 0.3; material.rho = 1000; % material
numScenarios = 1;

%% Construct FEA Solver
solver = feaClass(brep,numElements,material,vectorize,numScenarios, ...
    interpolation,penaltyStruct); % call superclass

solver = solver.fixEdge(6);
solver = solver.applyYForceOnEdge(3,-1e5);

solver = solver.preProcess(); % FEA pre-processing

%% Objective and Constraints
objective = densityComplianceElasticity(solver);

volumeFraction = 0.5;
constraints  = {volume(solver, volumeFraction)};

% manufacturing constraints
rmin = 1.5;
mfgConstraints = {
    minimumFeatureSize_dist(solver, rmin)
    symmetry_density(solver,1) % 0: x-dir, 1: y-dir
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
topopt.m_solver.plotDeformation();
topopt.m_solver.plotVonMisesStress();
topopt.m_solver.plotPrincipalStress();

%% Save Individual Figures
if exportImages
    saveAll(folder);%#ok
end

%% Export STL
if exportSTL
    thickness = 0.1;
    topopt.exportSTL(example_name, thickness);
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

