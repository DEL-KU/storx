clear; close all; format compact; format long

%% General parameters
vectorize = true;
exportImages = false;

%% File path
p = mfilename("fullpath");
[path,example_name,~] = fileparts(p);

disp("==================================");
disp(['Running ',example_name])

%% Problem definition
material.rho = 1;
material.mu = 1; % viscosity
numElements = 2000;

fem = fea2d_fluid('PipeBend.brep',numElements,material);

% inlet
Uin = 1;
fem = fem.fixUOfEdge(7,Uin);
fem = fem.fixVOfEdge(7,0);
% outlet
fem = fem.fixPOfEdge(2,0);
fem = fem.fixUOfEdge(2,0);

% no-slip top bottom
fem = fem.fixUOfEdge([1,3,5],0);
fem = fem.fixVOfEdge([1,3,5],0);
% no-slip left tight
fem = fem.fixUOfEdge([4,6,8],0);
fem = fem.fixVOfEdge([4,6,8],0);
%% Export
if exportImages
    % Make directory
    folder = [path '/../result/fea2d/example' '-' example_name '/']; %#ok
    mkdir(folder)
    cd(folder)
    delete 'log.txt'
    diary 'log.txt'
end

%% Solve
fem = fem.preProcess();
fem = fem.solve();
fem = fem.postProcess();

%% Output
fem.printFluidResults();

%% Plot
fem.plotGeometryWithLabels();
fem.plotMesh();
fem.plotBoundaryCondition();
fem.plotVelocity();
fem.plotPressure();

%% Save
if exportImages 
    saveAll(folder);%#ok
 end

%% Plot Combined Figures
ex_title = strjoin({'Fluid ','Example',example_name},' ');
combineFigures(ex_title);
if exportImages 
    saveAll(folder);%#ok
 end

cd(path)