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
material.k = 80;
numElements = 3200;
L = 0.2;
H = 0.1;
h = 0.07;
nFins = 5;
nVertices = 8 + 4*(nFins-2);
nEdges = nVertices;
t = L/(2*nFins-1);
v = [0 0; L 0; L H;L-t H;L-t H-h ];
xLoc = L-t;
for n = 1:nFins-2
    v(end+1,:) = [xLoc-t H-h];
    v(end+1,:) = [xLoc-t H];
    v(end+1,:) = [xLoc-2*t H];
    v(end+1,:) = [xLoc-2*t H-h];
    xLoc = xLoc - 2*t;
end
v(end+1,:) = [xLoc-t H-h];
v(end+1,:) = [xLoc-t H];
v(end+1,:) = [xLoc-2*t H];
heatFin.vertices = v';
heatFin.segments = zeros(4,nEdges);
heatFin.segments(1,:) = 1;
heatFin.segments(2,:) = 1:nEdges;
heatFin.segments(3,1:nEdges) = [2:nEdges 1];

fem = thermalClass(heatFin,numElements,material,vectorize);
fem = fem.fixEdge(2:nEdges,23);
fem = fem.applyFlux(1,1000);

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