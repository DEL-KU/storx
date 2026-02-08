clear; close all; format compact; format long
thermalClass = @fea2d_thermal;
elasticityClass = @fea2d_elasticity;
thermoelasticityClass = @fea2d_thermoelasticity;


%% General parameters
vectorize = true;
exportImages = false;

%% File path
p = mfilename("fullpath"); 
[path,example_name,~] = fileparts(p);

disp("==================================");
disp(['Running ',example_name])

%% Problem definition
brep = 'SquareSplit.brep';
numElements = 2000;
material.k = 1;
material.E = 2e11;  material.nu = 0.28; material.rho = 1;
material.alpha = 1.2e-5;
TReference = 0;

fem_thermal = thermalClass(brep,numElements,material,vectorize);
fem_thermal = fem_thermal.fixEdge(1:6,1);
fem_thermal = fem_thermal.applyInternalHeat(-0.1);

fem_elasticity = elasticityClass(brep,numElements,material,vectorize);
fem_elasticity = fem_elasticity.fixEdge([4,6]);
fem_elasticity = fem_elasticity.applyYForceOnEdge(2,1);

fem = thermoelasticityClass(fem_thermal,fem_elasticity, ...
    TReference,brep,numElements,material,vectorize);
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
fem.m_thermalSolver.printThermalResults();
fem.m_elasticitySolver.printElascticityResults();

%% Plot
fem.plotGeometryWithLabels();
fem.plotBoundaryCondition();
fem.plotMesh();
fem.m_thermalSolver.plotTemperature();
fem.m_elasticitySolver.plotDeformation();
fem.m_elasticitySolver.plotVonMisesStress();

%% Save
if exportImages 
    saveAll(folder);%#ok
 end

%% Plot Combined Figures
ex_title = strjoin({'Thermo-Elastic ','Example',example_name},' ');
combineFigures(ex_title);
if exportImages 
    saveAll(folder);%#ok
 end

cd(path)