function fem = verticalBarFlux(options)
arguments
    options.vectorize (1,1) logical = true
    options.exportImages (1,1) logical = false
    options.numElements (1,1) double {mustBeInteger,mustBePositive} = 1000
    options.material (1,1) struct = struct('k',1)
end

close all; format compact; format long
thermalClass = @fea2d_thermal;


%% General parameters
vectorize = options.vectorize;
exportImages = options.exportImages;

%% File path
p = mfilename("fullpath"); 
[path,example_name,~] = fileparts(p);

disp("==================================");
disp(['Running ',example_name])

%% Problem definition
material = options.material; % thermal conductivity
numElements = options.numElements;
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
end
