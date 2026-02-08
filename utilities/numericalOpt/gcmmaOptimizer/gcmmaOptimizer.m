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

classdef gcmmaOptimizer
    properties(GetAccess = 'public', SetAccess = 'public')
        m_n;
        m_m;
        m_outeriter;
        m_epsimin;
        m_raa0eps;
        m_raaeps;
        m_raa0;
        m_raa;
        m_ai;
        m_ci;
        m_di;
        m_a0;
        m_a;
        m_c;
        m_d;
        m_low;
        m_upp;
        m_fapp;
        m_xold1;
        m_xold2;
        m_f0app;
    end
    methods
        %%
        function obj = gcmmaOptimizer(n,m,ai,ci,di,x0)
            obj.m_n = n; % # variables
            obj.m_m = m; % # constraints
            obj.m_ai = ai;
            obj.m_ci = ci;
            obj.m_di = di;
            % DEFAULT VALUES
            obj.m_outeriter = 0;
            obj.m_epsimin = 1e-7;
            obj.m_raa0eps = 0.0000001;%1e-6;
            obj.m_raaeps = obj.m_raa0eps;
            obj.m_raa0 = 1e-4;
            obj.m_raa = 1e-4*ones(m,1);
            obj.m_a0 = 1;
            obj.m_a = obj.m_ai*ones(m,1);
            obj.m_c = obj.m_ci*ones(m,1);
            obj.m_d = obj.m_di*ones(m,1);
            obj.m_low = zeros(n,1);
            obj.m_upp = zeros(n,1);
            obj.m_f0app = 0;
            obj.m_fapp = zeros(m,1);
            obj.m_xold1 = x0;
            obj.m_xold2 = x0;
        end
        %%
        function [obj,xmma] = outerUpdate(obj,xval,fx,dfdx,gx,dgdx,xmin,xmax)

            obj.m_outeriter = obj.m_outeriter + 1;
            [obj.m_low,obj.m_upp,obj.m_raa0,obj.m_raa] = ...
                asymp(obj.m_outeriter,obj.m_n,xval,obj.m_xold1,obj.m_xold2,xmin,xmax,obj.m_low,obj.m_upp, ...
                obj.m_raa0,obj.m_raa,obj.m_raa0eps,obj.m_raaeps,dfdx,dgdx);
            % Generate the subproblem
            [xmma,~,~,~,~,~,~,~,~,obj.m_f0app,obj.m_fapp] = ...
                gcmmasub(obj.m_m,obj.m_n,obj.m_outeriter,obj.m_epsimin,xval,xmin,xmax,obj.m_low,obj.m_upp, ...
                obj.m_raa0,obj.m_raa,fx,dfdx,gx,dgdx,obj.m_a0,obj.m_a,obj.m_c,obj.m_d);
        end
        %%
        function [obj,xmma] = innerUpdate(obj,xnew,fxnew,gxnew, xval,fx, dfdx, gx, dgdx, xmin,xmax)
            %%%% New values on the parameters raa0 and raa are calculated:
            [obj.m_raa0,obj.m_raa] = ...
                raaupdate(xnew,xval,xmin,xmax,obj.m_low,obj.m_upp,fxnew,gxnew, ...
                obj.m_f0app,obj.m_fapp,obj.m_raa0,obj.m_raa,obj.m_raa0eps, ...
                obj.m_raaeps,obj.m_epsimin);
            %%%% The GCMMA subproblem is solved with these new raa0 and raa:
            [xmma,~,~,~,~,~,~,~,~,obj.m_f0app,obj.m_fapp] = ...
                gcmmasub(obj.m_m,obj.m_n,obj.m_outeriter,obj.m_epsimin,xval,xmin,xmax,obj.m_low,obj.m_upp, ...
                obj.m_raa0,obj.m_raa,fx,dfdx,gx,dgdx,obj.m_a0,obj.m_a,obj.m_c,obj.m_d);
        end
        %%
        function state = conCheck(obj,fx,gx)
            state = concheck(obj.m_m,obj.m_epsimin,obj.m_f0app,fx,obj.m_fapp,gx);
        end

        %%
    end
end
