clear; close all; format compact; format long
%% General Parameters
vectorize = true;
uniformGrid = 1;
exportImages = false;
exportGIF = false;
exportSTL = false;
%% Solvers
feaClass = @fea2d_elasticity;
topoptClass = @density2d_elasticity;

%% File path
p = mfilename("fullpath");
[path,example_name,~] = fileparts(p);

disp("==================================");
disp(['Running ',example_name])

%% Optimizer Parameters
interpolation = 'simp';
update = 'OC';
maxNumIters = 300;
penaltyStruct = struct('min',3,'max',3,'inc',0.0);

%% Problem definition
brep = 'LBracketNoFilletMidLoad.brep'; % geometry
numElements = 3200; % mesh
material.E = 100e9; material.nu = 0.3; material.rho = 1000; % material
numScenarios = 1;
%% Construct FEA Solver
solver = feaClass(brep,numElements,material,vectorize,numScenarios, ...
    interpolation,penaltyStruct,uniformGrid); % call superclass

solver = solver.fixEdge([7,8,9]);
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
    };

%% Construct Optimizer
topopt = topoptClass(solver, ...
    objective,constraints,mfgConstraints, ...
    update, ...
    maxNumIters,exportGIF);

%% Make Directory
if exportImages || exportSTL || exportGIF
    folder = [path '/../result/example' '-' example_name '/']; %#ok
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
cd(path)