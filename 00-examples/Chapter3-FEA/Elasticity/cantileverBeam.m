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
brep = 'CantileverBeam.brep'; % geometry
numElements = 3200; % mesh
material.E = 100e9; material.nu = 0.3; material.rho = 1000; % material

% construct fea solver
fem = elasticityClass(brep,numElements,material,vectorize); % call superclass

fem = fem.fixEdge(5);
fem = fem.applyYForceOnEdge(2,-1e5);


%% Export
if exportImages
    % Make directory
    folder = [path '/../result/example' '-' example_name '/']; %#ok
    mkdir(folder)
    cd(folder)
    delete 'log.txt'
    diary 'log.txt'
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
cd(path)