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

classdef ocOptimizer < handle
    properties(GetAccess = 'public', SetAccess = 'public')
        m_lambda1; % lower Lagrange multiplier (e.g., 0)
        m_lambda2; % upper Lagrange multiplier (e.g., 1000)
        m_move; % move limit (e.g., 0.2)
        m_eta; % fine-tuning parameter
        m_xmin; % minimum value of solution x (e.g., 0.001)
        m_xmax; % maximum value of solution x (e.g., 1)
        m_gmax; % right-hand side of constraint (sum(x) <= m_gmax)
    end
    methods
        %%
        function obj = ocOptimizer(l1,l2,move,xmin,xmax,gmax,eta)
            obj.m_lambda1 = l1;
            obj.m_lambda2 = l2;
            obj.m_move = move;
            obj.m_xmin = xmin;
            obj.m_xmax = xmax;
            obj.m_gmax = gmax;
            if (nargin < 7),obj.m_eta = 0.5;else,obj.m_eta = eta;end
        end
        %%
        function [obj,xnew] = update(obj,x,dfdx)
            l1 = obj.m_lambda1;  l2 = obj.m_lambda2; move = obj.m_move;
            while (l2-l1)/(l1+l2) > 1e-3
                lmid = 0.5*(l2+l1);
                B = -dfdx./lmid;
                xnew = min(x+move,x.*(B.^obj.m_eta)); % move along gradient
                xnew = max(x-move,min(obj.m_xmax,xnew)); % ensure it is less than max.
                xnew = max(obj.m_xmin, xnew);% ensure it is larger than min.
                if sum(xnew,'all')-obj.m_gmax > 0
                    l1 = lmid;
                else
                    l2 = lmid;
                end
            end
        end

    end
end
