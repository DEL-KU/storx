function solver = gridFEA_gripper(options)
arguments
    options.vectorize (1,1) logical = true
    options.exportImages (1,1) logical = false
    options.numElements (1,1) double {mustBeInteger,mustBePositive} = 4000
    options.material (1,1) struct = struct('E',2e9,'nu',0.35,'rho',1300)
    options.force (1,1) double = 10
    options.numScenarios (1,1) double {mustBeInteger,mustBePositive} = 1
end
configureGraphics();

close all;format compact; format long
warning('off','all')

%% Solvers
feaClass = @fea2d_elasticity;

%% General Parameters
vectorize = options.vectorize;
exportImages = options.exportImages;

%% File Path
p = mfilename("fullpath");
[path,example_name,~] = fileparts(p);

disp("==================================");
disp(['Running ',example_name])

%% Problem Definition
brep = 'GripperComplex.brep'; % geometry
numElements = options.numElements; % mesh
material = options.material;
force = options.force; % N
numScenarios = options.numScenarios;
%% Construct FEA Solver
solver = feaClass(brep,numElements,material,vectorize,numScenarios);

solver = solver.fixEdge([5,6,11,12]);
solver = solver.applyXForceOnEdge(18,force);

solver = solver.preProcess(); % FEA pre-processing
solver = solver.solve();
solver = solver.postProcess();
%% Make Directory
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

%% Output
solver.printElascticityResults();
solver.plotGeometryWithLabels();
solver.plotMesh();
solver.plotBoundaryCondition();
solver.plotDeformation('faceted');
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
if exportImages
    diary off
end

cd(path)
end
