function fem = HeatFin(options)
arguments
    options.vectorize (1,1) logical = true
    options.exportImages (1,1) logical = false
    options.combineFigure (1,1) logical = false
    options.numElements (1,1) double {mustBeInteger,mustBePositive} = 3200
    options.material (1,1) struct = struct('k',80)
    options.length (1,1) double {mustBePositive} = 0.2
    options.baseHeight (1,1) double {mustBePositive} = 0.1
    options.finHeight (1,1) double {mustBePositive} = 0.07
    options.nFins (1,1) double {mustBeInteger,mustBePositive} = 5
end

configureGraphics();
close all; format compact; format long
thermalClass = @fea2d_thermal;


%% General parameters
vectorize = options.vectorize;
exportImages = options.exportImages;
combineFigure = options.combineFigure;

%% File path
p = mfilename("fullpath"); 
[path,example_name,~] = fileparts(p);

disp("==================================");
disp(['Running ',example_name])

%% Problem definition
material = options.material;
numElements = options.numElements;
L = options.length;
H = options.baseHeight;
h = options.finHeight;
nFins = options.nFins;
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
if combineFigure
    ex_title = strjoin({'Thermal ','Example',example_name},' ');
    combineFigures(ex_title);
    if exportImages
        saveAll(folder);%#ok
    end
end
if exportImages
    diary off
end

cd(path)
end
