clc;clear;  close all;format compact; format long
warning('off','all')

%% Solvers
feaClass = @fea2d_elasticity;
topoptClass = @density2d_elasticity;

%% General Parameters
vectorize = true;
exportImages = true;
exportGif = false;
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

%% Objective and Constraints
objective = densityComplianceElasticity(solver);

volumeFraction = 0.65;
constraints  = {volume(solver, volumeFraction)};

% manufacturing constraints
rmin = 1.5;
mfgConstraints = {
    minimumFeatureSize_dist(solver, rmin)
    retain_density(solver,[5,6,11,12])
    }; 
%% Construct Optimizer
topopt = topoptClass(solver, ...
    objective,constraints,mfgConstraints, ...
    update, ...
    maxNumIters,exportGif);

%% Make Directory
if exportImages
    folder = [path '/result/example_noRetain' '-' example_name '/']; %#ok
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
topopt.plotIsoSurface('Contour');
topopt.m_solver.plotDeformation();
topopt.m_solver.plotVonMisesStress();
topopt.m_solver.plotPrincipalStress();

%% Export STL
if exportSTL
    thickness = 10;
    topopt.exportSTL(example_name, thickness);
end

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
if exportImages
    diary off
end

cd(path)