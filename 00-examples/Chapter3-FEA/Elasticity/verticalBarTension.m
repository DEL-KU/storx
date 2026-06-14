clear; close all; format compact; format long
elasticityClass = @fea2d_elasticity;

%% General parameters
vectorize = true;
exportImages = false;

%% File path
p = mfilename("fullpath"); 
[path,example_name,~] = fileparts(p);

disp("==================================");
disp(['Running ',example_name])

%% Problem definition
material.E = 2e11;  material.nu = 0.28;
material.rho = 1;
numElements = 1000;
verticalBar.vertices = [0 0; 1 0; 1 10; 0 10]';
verticalBar.segments = [1 1 2 0 ;1 2 3 0;1 3 4 0;1 4 1 0]';
fem = fea2d_elasticity(verticalBar,numElements,material,vectorize);
fem = fem.fixYOfEdge(1);
fem = fem.applyYForceOnEdge(3,100);
disp('Exact answer for tensile load: delta = 5e-9, sigma = 100')
example_name =  'verticalBarTension';

%% Export
if exportImages
    % Make directory
    folder = [path '/../result/example' '-' example_name '/']; %#ok
    mkdir(folder)
    cd(folder)
    diary off
    logFile = fullfile(folder, 'log.txt');
    if exist(logFile, 'file')
        delete(logFile)
    end
    diary(logFile)
end

%% Solve
fem = fem.preProcess();
fem = fem.solve();
fem = fem.postProcess(true);

%% Output
fem.printElascticityResults();
fem.plotGeometryWithLabels();
fem.plotMesh();
fem.plotBoundaryCondition();
fem.plotDeformation();
fem.plotVonMisesStress();
fem.plotPrincipalStress();

%% Save
if exportImages 
    saveAll(folder);%#ok
 end

%% Plot Combined Figures
ex_title = strjoin({'Elasticity ','Example',example_name},' ');
combineFigures(ex_title);
if exportImages 
    saveAll(folder);%#ok
 end
if exportImages
    diary off
end

cd(path)