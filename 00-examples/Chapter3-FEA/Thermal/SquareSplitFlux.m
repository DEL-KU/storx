clear; close all; format compact; format long
thermalClass = @fea2d_thermal;


%% General parameters
vectorize = true;
exportImages = false;

%% File path
p = mfilename("fullpath"); 
[path,example_name,~] = fileparts(p);

disp("==================================");
disp(['Running ',example_name])

%% Problem definition
material.k = 1;
numElements = 10000;
fem = thermalClass('SquareSplit.brep',numElements,material,vectorize);
fem = fem.fixEdge(2,0);
fem = fem.applyFlux(5,1);      

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
fem = fem.postProcess();

%% Output
fem.printThermalResults();

%% Plot
fem.plotGeometryWithLabels();
fem.plotMesh();
fem.plotBoundaryCondition();
fem.plotTemperature();
%% Save
if exportImages 
    saveAll(folder);%#ok
 end

%% Plot Combined Figures
ex_title = strjoin({'Thermal ','Example',example_name},' ');
combineFigures(ex_title);
if exportImages 
    saveAll(folder);%#ok
 end
if exportImages
    diary off
end

cd(path)