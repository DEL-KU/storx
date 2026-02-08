clear; close all; format compact; format long
elasticityClass = @triFEA2d_elasticity;
%% General parameters
exportImages = false;

%% File path
p = mfilename("fullpath"); 
[path,example_name,~] = fileparts(p);

disp("==================================");
disp(['Running ',example_name])

%% Problem definition
brep = 'SquareWithHole.brep'; % geometry
numElements = 1000; % mesh
material.E = 2e9;  material.nu = 0.28; material.rho = 1;
% construct fea solver
fem = elasticityClass(brep,numElements,material); % call superclass

fem = fem.fixEdge(9);   
fem = fem.applyXForceOnEdge(7,1);


%% Export
if exportImages
    % Make directory
    folder = [path '/../result/example' '-' example_name '/']; %#ok
    mkdir(folder)
    cd(folder)
    delete 'log.txt'
    diary 'log.txt'
end

% Solve
fem = fem.preProcess();
fem = fem.solve();
fem = fem.postProcess(true);

% Output
fem.printElascticityResults();
fem.plotGeometryWithLabels();
fem.plotMesh();
fem.plotBoundaryCondition();
fem.plotDeformation();
fem.plotVonMisesStress();
% fem.plotPrincipalStress();

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