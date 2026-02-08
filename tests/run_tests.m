clc
%% REPRESENTATION
runtests('test_brep')
runtests('test_gridMesher')

%% FEA 2D
runtests('test_fea2d_elasticity')
runtests('test_fea2d_thermal')
runtests('test_fea2d_thermoelasticity')

%% PARAMETRIC OPTIMIZATION
runtests('test_parametricOptimization')

%% TOPOLOGY OPTIMIZATION
% Density
runtests('test_topopt2d_density')

% Levelset
runtests('test_topopt2d_levelset')

% % Evolutionary
runtests('test_topopt2d_evolutionary')

% % PareTO
runtests('test_topopt2d_pareto')