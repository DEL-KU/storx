% This function is a modification of 'reinitTest.m' code in 'ToolboxLS':
%   #######################################################################
%   Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
%   This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%   Ian Mitchell, 2/14/04
%   #######################################################################
% Amir M. Mirzendehdel, 8/24/2022

function sdf = reinit2D(lsf,dx,boundingBox,accuracy)
if (nargin < 4), accuracy = 'medium'; end
%--------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
addpath(genpath('.\ToolboxLS'));
%--------------------------------------------------------------------------
% Integration parameters.
tMax = 1.0;                  % End time.
tStep = 0.1;                 % Period at which plot should be produced.
t0 = 0;                      % Start time.
small = 100 * eps;           % Convergence relative to tMax
%--------------------------------------------------------------------------
% Create the grid.
g.dim = 2;
g.min=boundingBox(:,1)-0.5*dx;
g.max=boundingBox(:,2)+0.5*dx;
g.dx = dx;
g.bdry = @addGhostExtrapolate;
g = processGrid(g);
%--------------------------------------------------------------------------
% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.5, 'stats', 'off');
%--------------------------------------------------------------------------
% Choose approximations at appropriate level of accuracy.
switch(accuracy)
 case 'low'
  derivFunc = @upwindFirstFirst;
  integratorFunc = @odeCFL1;
 case 'medium'
  derivFunc = @upwindFirstENO2;
  integratorFunc = @odeCFL2;
 case 'high'
  derivFunc = @upwindFirstENO3;
  integratorFunc = @odeCFL3;
 case 'veryHigh'
  derivFunc = @upwindFirstWENO5;
  integratorFunc = @odeCFL3;
 otherwise
  error('Unknown accuracy level %s', accuracy);
end
%--------------------------------------------------------------------------
data = lsf';
% Set up spatial approximation scheme.
schemeFunc = @termReinit;
schemeData.grid = g;
schemeData.derivFunc = derivFunc;
schemeData.initial = data;
% Use the subcell fix by default.
schemeData.subcell_fix_order = 1;
%--------------------------------------------------------------------------
% Loop until tMax (subject to a little roundoff).
tNow = t0;
while(tMax - tNow > small * tMax)
  % Reshape data array into column vector for ode solver call.
  y0 = data(:);
  % How far to step?
  tSpan = [ tNow, min(tMax, tNow + tStep) ];
  % Take a timestep.
   [t,y] = integratorFunc( schemeFunc, tSpan, y0, ...
       integratorOptions, schemeData);
  tNow = t(end);
  % Get back the correctly shaped data array
  data = reshape(y, g.shape);
end
%--------------------------------------------------------------------------
% Return SDF 
sdf = data';
