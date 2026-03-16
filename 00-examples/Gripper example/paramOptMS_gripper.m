clc;clear;  close all;format compact; format long
warning('off','all')

%% General Parameters
exportImages = true;
exportGif = false;

%% File Path
p = mfilename("fullpath"); 
[path,example_name,~] = fileparts(p);

%% Export
if exportImages
    % Make directory
    folder = [path '/result/example' '-' example_name '/']; %#ok
    mkdir(folder)
    cd(folder)
    delete 'log.txt'
    diary 'log.txt'
end

disp("==================================");
disp(['Running ',example_name])

%% Problem Definition
params0.value = [20,-15, 4.0,10,-60,5.0];
params0.lb = [10,  -10,   1.0, 5,  -40,   1.0 ];
params0.ub = [ 25,  -25,   10.0,10,  -55,  4.0 ];

objective = 'compliance'; % objective

constraints.area = 3.1e3; % constraint value
constraints.type = 'ineq'; % constraint type: 'eq' or 'ineq'
%% Construct Optimizer
brepHandle = @createGeom;
solverHandle = @createProblem;
terminationTolerance = 1e-6;
finiteDifferenceStepSize = 1e-6;

% Optimization method:
%  - RS: Random Search
%  - FD: Finite Difference
%  - MS: Multi-Start
%  - GS: Global Search
method = "RS";

parOpt = parameterOpt2d(brepHandle,solverHandle,params0, ...
    objective,constraints, ...
    terminationTolerance,finiteDifferenceStepSize,method,exportGif);

parOpt = parOpt.setNumberOfRandomSearchSamples(10, 100);
%% Optimize
parOpt = parOpt.optimize();

%% Save
if exportImages 
    saveAll(folder);%#ok
 end

%% Plot Combined Figures
ex_title = strjoin({'Parametric Shape Opt. ','Example',example_name},' ');
combineFigures(ex_title);
if exportImages 
    saveAll(folder);%#ok
 end
cd(path)

%% Create Problem
function solver = createProblem(brep)
vectorize = true;
numElements = 10000; % mesh
material.E = 2e9; material.nu = 0.35; material.rho = 1300; % material
force = 10; % N
numScenarios = 1;
solver = fea2d_elasticity(brep,numElements,material,vectorize,numScenarios); % call superclass

solver = solver.fixEdge([5,6,11,12]);
solver = solver.applyXForceOnEdge(18,force);

end

%% Create Geometry from Parameters
function geom = createGeom(params)
%CREATEGEOM  Parametric B-rep with two additional internal circular holes.
%
% params = [x3, y3, r3,  x4, y4, r4]
%   hole #3 center = (x3,y3), radius = r3
%   hole #4 center = (x4,y4), radius = r4
%
% Notes:
% - Keeps your original outer boundary + original internal holes unchanged.
% - Adds two more internal holes using your (-1) virtual edges + two (2) arcs.
% - Uses diametrically-opposed points to represent each full circle via two arcs.
%
% You should enforce bounds in your optimizer; this function assumes params are valid.

arguments
    params (1,6) double
end

x3 = params(1);  y3 = params(2);  r3 = params(3);
x4 = params(4);  y4 = params(5);  r4 = params(6);

% --- Base vertices (IDs 1..16) ---
V = [ ...
   -16     0.96   40.3   24.09   14.56    3.66   -7.67     0    -0.8   -16    19     0    21.5   16.5     0      0 ; ...
  -100  -100    -17.68   22.08   16.75    7.11   -2.27  -28.22  -82    -82    20     0     20     20      4     -4 ] ;

% --- Base segments (IDs reference above vertices) ---
S = [ ...
    1   1   2    0
    1   2   3    0
    1   3   4    0
   -1   4  13    0
    2  13  14   11
    2  14  13   11
   -1  13   4    0
    2   4   5  -11
    1   5   6    0
   -1   6  15    0
    2  15  16   12
    2  16  15   12
   -1  15   6    0
    2   6   7  -12
    1   7   8    0
    1   8   9    0
    1   9  10    0
    1  10   1    0 ] ;

% --- New hole centers (parametric) ---
C3 = [x3; y3];
C4 = [x4; y4];

% Define two opposite points on each circle (along +x / -x from center)
P3a = C3 + [ r3; 0];
P3b = C3 + [-r3; 0];
P4a = C4 + [ r4; 0];
P4b = C4 + [-r4; 0];

% Append vertices: for each new hole add (center, pointA, pointB)
% New IDs:
%   hole #3: c3=17, a3=18, b3=19
%   hole #4: c4=20, a4=21, b4=22
V = [V, C3, P3a, P3b, C4, P4a, P4b];
c3 = 17; a3 = 18; b3 = 19;
c4 = 20; a4 = 21; b4 = 22;

% Virtual anchors (existing boundary vertices) to "register" holes to the face
anchor3 = 3;   % near top-right-ish
anchor4 = 8;   % mid-left-ish
% (You can change anchors; they do not affect hole location, only topology encoding.)

S_add3 = [ ...
   -1  anchor3  a3   0
    2  a3       b3   c3
    2  b3       a3   c3
   -1  a3  anchor3   0 ];

S_add4 = [ ...
   -1  anchor4  a4   0
    2  a4       b4   c4
    2  b4       a4   c4
   -1  a4  anchor4   0 ];

S = [S; S_add3; S_add4];

geom.vertices = V;      % 2-by-N
geom.segments = S.';    % 4-by-M
end
