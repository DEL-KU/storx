clear; close all; format compact; format long
elasticityClass = @triFEA2d_elasticity;
%% General parameters
exportImages = true;

%% File path
p = mfilename("fullpath");
[path,example_name,~] = fileparts(p);

disp("==================================");
disp(['Running ',example_name])

%% Problem definition
brep = 'GripperComplex.brep'; % geometry
numElements = 6000; % mesh
material.E = 2e9; material.nu = 0.35; material.rho = 1300; % material
force = 10; % N
numScenarios = 1;
class = 'PlaneStress';
order = 'Linear';

% construct fea solver
fem = elasticityClass(brep,numElements,material,class,order); % call superclass

fem = fem.fixEdge([5,6,11,12]);
fem = fem.applyXForceOnEdge(18,force);

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

%% Save
if exportImages
    saveAll(folder);%#ok
end

%% Plot Combined Figures
% ex_title = strjoin({'Elasticity ','Example',example_name},' ');
% combineFigures(ex_title);
% if exportImages
%     saveAll(folder);%#ok
% end
% cd(path)