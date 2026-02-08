%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Matlab code was written by:                                          %
% - Amir M. Mirzendehdel, Aerospace Engineering Department, KU              %
% - Krishnan Suresh, Mechanical Engineering Department, UW-Madison          %
%                                                                           %
% Please send your comments to: amirzend@ku.edu                             %
%                                                                           %
% The code is intended for educational purposes and theoretical details     %
% are discussed in the textbook:                                            %
% Introduction to Shape and Topology Optimization using MATLAB              %
%                                                                           %
% Disclaimer:                                                               %
% The authors reserves all rights but do not guaranty that the code is      %
% free from errors. Furthermore, we shall not be liable in any event        %
% caused by the use of the program.                                         %
%                                                                           %
% License:                                                                  %
% This software is used, copied and distributed under the licensing         %
% agreement contained in the file LICENSE in the top directory of           %
% the distribution.                                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef PlotId
   properties
      %% boundary represenation
      brep = 1;
      polygon = 100;
      brep_optimized = 500;
      %% mesh represenation
      mesh = 1000;
      %% boundary conditions
      loading = 2000;
      %% shape/topology optimization
      initial_holes = 3000;
      design = 4000;
      ls_dsdx = 4100;
      ls_velcity = 4200;
      lsf = 4300;
      isosurface_contour = 5000;
      isosurface_lsf = 6000;
      %% composite optimization
      composite_discrete= 7000;
      composite_streamline = 8000;
      %% elasticity
      deformation = 9000;
      von_mises = 10000;
      principal_stress = 11000;
      %% thermal
      temperature = 12000;
      %% fluid
      velocity = 13000;
      pressure = 14000;
      %% convergence
      convergence = 15000;
      pareto_front = 16000;
      %% geometries
      dxf_lines = 17000;
      geom_stl = 18000;
   end
end

