%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discription:                                                              %
% This is a class for functionals used in shape and               %
% topology optimization as objectives and constraints.                      %
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

classdef volume < functional
    properties (GetAccess = 'public', SetAccess = 'protected')

    end

    methods
        %% CONSTRUCTOR
        function obj = volume(solver, ub)

            % if upper bound is provided set the value, otherwise set it to NaN
            upper_bound = NaN;
            if (nargin > 1)
                upper_bound = ub;
            end

            % constructor based on superclass functional
            % solver can be any instance of simulation2d class
            obj = obj@functional(solver, upper_bound);
        end

        % evaluate the compliance volume at design and state variables queried from the elasticity solver
        % input: obj, ubound (optional)
        % output: obj, volume
        % ubound is optional, if not provided, the functional is considered as objective and returns fx
        % if ubound is provided, the functional is considered as constraint and returns gx := fx/ubound - 1 <= 0
        function [obj,value] = evaluate(obj)
            obj.m_fx = sum(obj.m_solver.m_design(:))/(obj.m_solver.m_numExistingElems);
            if (obj.m_numEvaluations == 0), obj.m_fx0 = obj.m_fx; end

            if (obj.m_isConstraint)
                value = obj.m_fx/obj.m_upperBound - 1; % as constraint
            else
                value = obj.m_fx; % as objective
            end
            obj.m_numEvaluations = obj.m_numEvaluations + 1; % update number of evaluations
        end

        % gradient of the compliance at design and state variables queried from the solver.
        % additional adjoint problems are solved in this function to obtain the lagrnage multipliers
        % input: obj, ubound (optional)
        % output: obj, functional gradient w.r.t. design variables
        % ubound is optional, if not provided, the functional is considered as objective and returns df/dx
        % if ubound is provided, the functional is considered as constraint and returns dgx/dx := dfx/dx/ubound
        function [obj,grad] = gradient(obj)
            obj.m_dfdx = 1/(obj.m_solver.m_numExistingElems) * ones(size(obj.m_solver.m_design));

            if (obj.m_isConstraint)
                grad = obj.m_dfdx/obj.m_upperBound;
            else
                grad = obj.m_dfdx;
            end

        end
    end
end
