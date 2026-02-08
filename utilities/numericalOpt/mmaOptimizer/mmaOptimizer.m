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

classdef mmaOptimizer
    properties(GetAccess = 'public', SetAccess = 'public')
        m_n; % # variables
        m_m; % # constraints
        m_iter;
        m_xmamieps; m_epsimin;  m_raa0; m_move; m_albefa;
        m_asyminit; m_asymdec;  m_asyminc;
        m_ai;   m_ci;   m_di;
        m_a;    m_c;    m_d;    m_y;    m_b;
        m_lam;  m_mu;   m_s;
        m_low;  m_upp;
        m_alpha;    m_beta;
        m_p0;   m_q0;
        m_pij;  m_qij;
        m_grad; m_hess;
        m_xold1;    m_xold2;
        m_z;
    end
    methods
        %%
        function obj = mmaOptimizer(n,m,ai,ci,di,x0)
            obj.m_n = n; % # variables
            obj.m_m = m; % # constraints
            obj.m_ai = ai;
            obj.m_ci = ci;
            obj.m_di = di;
            % DEFAULT VALUES
            obj.m_iter = 0;
            obj.m_xmamieps = 1e-5;
            obj.m_epsimin = 1e-7;
            obj.m_raa0 = 1e-5;
            obj.m_move = 0.2;
            obj.m_albefa = 0.1;
            obj.m_asyminit = 0.5;
            obj.m_asymdec = 0.7;
            obj.m_asyminc = 1.2;
            obj.m_a = obj.m_ai*ones(m,1);
            obj.m_c = obj.m_ci*ones(m,1);
            obj.m_d = obj.m_di*ones(m,1);
            obj.m_y = zeros(m,1);
            obj.m_lam = zeros(m,1);
            obj.m_mu = zeros(m,1);
            obj.m_s = zeros(m,1);
            obj.m_low = zeros(n,1);
            obj.m_upp = zeros(n,1);
            obj.m_alpha = zeros(n,1);
            obj.m_beta = zeros(n,1);
            obj.m_p0 = zeros(n,1);
            obj.m_q0 = zeros(n,1);
            obj.m_pij = zeros(m*n,1);
            obj.m_qij = zeros(m*n,1);
            obj.m_b = zeros(m,1);
            obj.m_grad = zeros(m,1);
            obj.m_hess = zeros(m*m,1);
            obj.m_xold1 = x0;
            obj.m_xold2 = x0;
        end
        %%
        function [obj,xmma] = update(obj,x,fx,dfdx,gx,dgdx,xmin,xmax)
            [xmma,~,~,~,~,~,~,~,~,obj.m_low,obj.m_upp] = ...
                mmasub(obj.m_m,obj.m_n,obj.m_iter,x,xmin,xmax,obj.m_xold1,obj.m_xold2, ...
                fx,dfdx,gx,dgdx,obj.m_low,obj.m_upp, ...
                obj.m_ai, obj.m_a, obj.m_c, obj.m_d);
            obj.m_xold2 = obj.m_xold1;
            obj.m_xold1 = x;
            obj.m_iter  = obj.m_iter + 1;
        end
    end
end
