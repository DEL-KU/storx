clc;clear;  close all;format compact; format long
warning('off','all')

%% Solvers
feaClass = @fea2d_elasticity;
topoptClass = @standardHJ2d_elasticity;

%% General Parameters
uniformGrid = 1; % needed for the Hamilton-Jacobi solver
vectorize = true;
exportImages = true;
exportGif = true;
exportSTL = true;

%% File Path
p = mfilename("fullpath");
[path,example_name,~] = fileparts(p);

disp("==================================");
disp(['Running ',example_name])

%% Optimizer Parameters
interpolation = 'none';
maxNumIters = 600;
penaltyStruct = struct('min',1,'max',1,'inc',0);

%% Problem Definition
brep = 'GripperComplex.brep'; % geometry
numElements = 10000; % mesh
material.E = 2e9; material.nu = 0.35; material.rho = 1300; % material
force = 10; % N
numScenarios = 1;
%% Construct FEA Solver
solver = feaClass(brep,numElements,material,vectorize,numScenarios, ...
    interpolation,penaltyStruct,uniformGrid); % call superclass

solver = solver.fixEdge([5,6,11,12]);
solver = solver.applyXForceOnEdge(18,force);

solver = solver.preProcess(); % FEA pre-processing

%% Objective and Constraints
objective = standardHJComplianceElasticity(solver);

volumeFraction = 0.65;
constraints  = {volume(solver, volumeFraction)};

% manufacturing constraints
mfgConstraints = {
    minimumFeatureSize_conv(solver)
    retain_levelset(solver,[5,6,11,12,18])
    }; 
%% Construct Optimizer
nHolesX = 2; nHolesY = 4; r0 = 0.5;
topopt = topoptClass(solver, ...
                objective,constraints,mfgConstraints, ...
                nHolesX,nHolesY,r0, ...
                maxNumIters,exportGif);
%% Make Directory
if exportImages
    folder = [path '/result/example' '-' example_name '/']; %#ok
    name = ['numElem' num2str(numElements) '-' 'vf' num2str(volumeFraction)];
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
ex_title = strjoin({'Level-set SO for Elasticity ','Example',example_name},' ');
combineFigures(ex_title);
if exportImages
    saveAll(folder);%#ok
end
cd(path)