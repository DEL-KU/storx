clc;clear;  close all;format compact; format long
warning('off','all')

%% Solvers
feaClass = @fea2d_elasticity;
topoptClass = @density2d_elasticity;

%% General Parameters
vectorize = true;
exportImages = true;
exportGif = false;
exportSTL = true;

%% File Path
p = mfilename("fullpath");
[path,example_name,~] = fileparts(p);

disp("==================================");
disp(['Running ',example_name])

%% Optimizer Parameters
interpolation = 'simp';
update = 'MMA';
maxNumIters = 300;
penaltyStruct = struct('min',3,'max',5,'inc',0.01);

%% Problem Definition
brep = 'GripperComplex.brep'; % geometry
numElements = 80000; % mesh
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


% local volume constraint
volumeFraction = 0.65;
localRadius = 10;
localPNorm = 16;
constraints = {localVolume(solver, localRadius, localPNorm, volumeFraction)};

% manufacturing constraints
% min. feat. size filter
rmin = 2.5;

% Heaviside projection
mfgConstraints = {
    minimumFeatureSize_dist(solver, rmin)
    physicalDensity(solver)
    retain_density(solver,[5,6,11,12])
    }; 
%% Construct Optimizer
topopt = topoptClass(solver, ...
    objective,constraints,mfgConstraints, ...
    update, ...
    maxNumIters,exportGif);

%% Make Directory
if exportImages
    folder = [path '/result/example_localVol' '-' example_name '/']; %#ok
    name = [update '-' 'numElem' num2str(numElements) '-' 'vf' num2str(volumeFraction)];
    folder = [folder name '/'];
    mkdir(folder)
    cd(folder)
    delete 'log.txt'
    diary 'log.txt'
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
% ex_title = strjoin({example_name,'Combined '},' ');
% combineFigures(ex_title);
% if exportImages
%     saveAll(folder);%#ok
% end
% cd(path)