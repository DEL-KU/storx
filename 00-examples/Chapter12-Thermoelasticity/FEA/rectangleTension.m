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
material.E = 2e11;  material.nu = 0; material.rho = 1000;
material.k = 45;
material.alpha = 1e-5;
TReference = 0;

numElements = 100;
xmax = 3;
domain.vertices = [0 0; xmax 0; xmax 1; 0 1]';
domain.segments = [1 1 2 0 ;1 2 3 0;1 3 4 0;1 4 1 0]';

fem_thermal = thermalClass(domain,numElements,material,vectorize);
fem_thermal = fem_thermal.fixEdge(4,0); % fixed T
Q = 1000;
fem_thermal = fem_thermal.applyFlux(2,Q); % flux

tau = 0;
fem_elasticity = fea2d_elasticity(domain,numElements,material,vectorize);
fem_elasticity = fem_elasticity.fixXOfEdge(4);
fem_elasticity = fem_elasticity.fixYOfEdge([1,3]);
fem = thermoelasticityClass(fem_thermal,fem_elasticity, ...
    TReference,domain,numElements,material,vectorize);

% Analytical solution:
E = material.E; nu = material.nu; k = material.k;alpha = material.alpha;
lambda = (E*nu)/((1+nu)*(1-2*nu));
mu = E/(2*(1+nu));

u_max = Q/k*alpha * (lambda+mu)/(lambda+2*mu) * xmax^2 + tau/(lambda+2*mu)*xmax;
T_max = Q/k*xmax;
disp(['Exact answer for tensile load: u_max = ' num2str(u_max) ',T_max = ' num2str(T_max)])

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
ex_title = strjoin({'Thermo-Elastic:','Example',example_name},' ');
combineFigures(ex_title);
if exportImages 
    saveAll(folder);%#ok
 end

cd(path)