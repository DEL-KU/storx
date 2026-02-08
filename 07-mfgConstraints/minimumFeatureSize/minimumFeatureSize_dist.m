%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discription:                                                              %
% This is a class for imposing minimum feature size constraints             %
% based on distance between elements used in topology optimization          % 
% (primarily used in density-based TO such as SIMP).                        %
%                                                                           %
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

classdef  minimumFeatureSize_dist < mfgConstraints
    properties (GetAccess = 'public', SetAccess = 'protected')
        m_rmin; % minimum feature size
        m_H; % Heaviside filter
        m_Hs; % sum of Heaviside filter
    end
    methods
        %% CONSTRUCTOR
        function obj = minimumFeatureSize_dist(solver,rmin)
            % check if solver is valid
            if (~isa(solver, 'simulation2d')), error('solver must be an instance of simulation2d class!');end % check if solver is valid

            % constructor based on superclass
            obj = obj@mfgConstraints(solver);
            obj.m_rmin = rmin; % assign minimum feature size
            obj = obj.setupMinimumFeatureSizeFilter(rmin);
        end


        function [filteredDesign] = filterDesign(obj, design)
            % filter the design variables
            % input: obj, design variables
            % output: obj, filtered design variables
            filteredDesign(:) = obj.m_H*(obj.m_solver.m_existingElems(:).*design(:))./obj.m_Hs;
            filteredDesign = reshape(filteredDesign, size(design));
            filteredDesign = filteredDesign .* obj.m_solver.m_existingElems;
        end

        function [filteredSensitivity] = filterSensitivity(obj, design, sensField)
            % filter the sensitivity fields
            % input: obj, sensitivity fields
            % output: obj, filtered sensitivity fields
            filteredSensitivity(:) = obj.m_H*(design(:).*sensField(:))./obj.m_Hs./max(1e-3,design(:));
            filteredSensitivity = reshape(filteredSensitivity, size(design));
            filteredSensitivity = filteredSensitivity .* obj.m_solver.m_existingElems;
        end

    end
    methods (Access = private)
        function obj = setupMinimumFeatureSizeFilter(obj,rmin)
            % setup minimum feature size filter
            % input: obj, minimum feature size
            % output: obj
            nx = obj.m_solver.m_nx;
            ny = obj.m_solver.m_ny;
            iH = ones(nx*ny*(2*(ceil(rmin)-1)+1)^2,1);
            jH = ones(size(iH));
            sH = zeros(size(iH));
            k = 0;
            for i1 = 1:nx
                for j1 = 1:ny
                    e1 = (i1-1)*ny+j1;
                    for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nx)
                        for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),ny)
                            e2 = (i2-1)*ny+j2;
                            k = k+1;
                            iH(k) = e1;
                            jH(k) = e2;
                            if (~obj.m_solver.m_existingElems(j2,i2))
                                val = 1e-9;
                            else
                                val = rmin-sqrt((i1-i2)^2+(j1-j2)^2);
                            end
                            sH(k) = max(0,val);
                        end
                    end
                end
            end
            obj.m_H = sparse(iH,jH,sH);
            obj.m_Hs = sum(obj.m_H,2);
        end
    end
end
