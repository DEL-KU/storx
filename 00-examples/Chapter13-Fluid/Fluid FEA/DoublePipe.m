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

fem = fea2d_fluid('DoublePipe.brep',numElements,material);

% inlet
Uin = 1;
fem = fem.fixUOfEdge([9,11],Uin);
fem = fem.fixVOfEdge([9,11],0);
% outlet
fem = fem.fixPOfEdge([3,5],0);
fem = fem.fixVOfEdge([3,5],0);

% no-slip top bottom
fem = fem.fixUOfEdge([1,7],0);
fem = fem.fixVOfEdge([1,7],0);
% no-slip left tight
fem = fem.fixUOfEdge([2,4,6,8,10,12],0);
fem = fem.fixVOfEdge([2,4,6,8,10,12],0);

%% Export
if exportImages
    % Make directory
    folder = [path '/../result/fea2d/example' '-' example_name '/']; %#ok
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

if exportImages
    diary off
end

cd(path)