%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                              %
% This class implements a 2D evolutionary-based topology optimization using %
% evolutionary optimization. It inherits from the topopt2d class and        %
% implements the methods necessary to solve the Bi-directional evolutionary %
% structural optimization (BESO) problems.                                  %
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

classdef (Abstract) beso2d  < eso2d
    methods (Abstract)
        obj = evaluate(obj)
        obj = gradient(obj)
        obj = saveHistory(obj)
    end
    methods
        function obj = beso2d(solver,objective,constraints, mfgConstraints, ...
                volDecrement,exportGIF,testMode)

            % construct
            obj = obj@eso2d(solver,objective,constraints, mfgConstraints, ...
                volDecrement,exportGIF,testMode); % call superclass

            % history
            obj.m_history.constraint.volFrac = zeros(obj.m_maxNumIters,1);
            obj.m_history.state.deformation = zeros(obj.m_maxNumIters,obj.m_solver.m_numScenarios);
            obj.m_history.state.vonMises = zeros(obj.m_maxNumIters,obj.m_solver.m_numScenarios);
        end

        %% UPDATE DESIGN
        function obj = update(obj, volFrac)
            % Minimal BESO update (OC-style thresholding via bisection)

            % sensitivities (collapse scenarios if present)
            dfdx = -obj.m_dfdx;                            % ny x nx x (nScen)
            if ndims(dfdx) == 3, dfdx = mean(dfdx,3); end
            dc = -dfdx;                                   % "benefit" for keeping solid
            dv = obj.m_solver.m_ve;                       % uniform elem volume
            phi = dc / max(dv, eps);                      % efficiency ratio

            % mask & targets
            mask   = logical(obj.m_solver.m_existingElems);
            nExist = nnz(mask);
            nTar   = max(0, min(nExist, round(volFrac*nExist)));

            % trivial cases
            if nExist == 0
                obj.m_x(:) = 0;
                obj.m_solver = obj.m_solver.setDesign(obj.m_x);
                return
            end

            % bisection on threshold th so that #solids == nTar
            pv   = phi(mask);                      % work only on selectable elems
            l1   = min(pv);
            l2   = max(pv);
            for it = 1:50
                th = 0.5*(l1+l2);
                x_new = false(size(phi));
                x_new(mask) = pv >= th;            % keep if phi >= th
                nSol = nnz(x_new(mask));
                if nSol > nTar
                    l1 = th;                       % too many solids -> raise th
                elseif nSol < nTar
                    l2 = th;                       % too few solids -> lower th
                else
                    break
                end
                if (l2-l1) <= 1e-12*max(1,abs(l2)), break; end
            end

            % finalize
            obj.m_x = double(x_new);
            obj.m_x(~mask) = 0;                    % enforce forbidden cells void
            obj.m_solver = obj.m_solver.setDesign(obj.m_x);

            obj.m_tau = obj.findContourValueWithVolumeFraction(volFrac);
        end
    end
end
