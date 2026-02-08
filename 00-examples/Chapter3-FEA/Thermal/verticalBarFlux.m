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
material.k = 1; % thermal conductivity
numElements = 1000;
verticalBar.vertices = [0 0; 1 0; 1 10; 0 10]';
verticalBar.segments = [1 1 2 0 ;1 2 3 0;1 3 4 0;1 4 1 0]';
fem = thermalClass(verticalBar,numElements,material,vectorize);
fem = fem.fixEdge(1,0);
fem = fem.applyFlux(3,1);
disp('Exact answer: TMax = 10')
example_name = 'verticalBar';
        

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
cd(path)