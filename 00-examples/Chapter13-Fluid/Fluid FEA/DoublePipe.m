function fem = DoublePipe(options)
arguments
    options.exportImages (1,1) logical = false
    options.brep (1,:) char = 'DoublePipe.brep'
    options.numElements (1,1) double {mustBeInteger,mustBePositive} = 2000
    options.material (1,1) struct = struct('rho',1,'mu',1)
    options.inletVelocity (1,1) double = 1
end

close all; format compact; format long

%% General parameters
exportImages = options.exportImages;

%% File path
p = mfilename("fullpath"); 
[path,example_name,~] = fileparts(p);

disp("==================================");
disp(['Running ',example_name])

%% Problem definition
material = options.material;
numElements = options.numElements;

fem = fea2d_fluid(options.brep,numElements,material);

% inlet
Uin = options.inletVelocity;
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
end
