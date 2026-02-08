%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discription:                                                              %
% This is an abstract class for functionals used in shape and               %
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

classdef (Abstract)  functional
    properties (GetAccess = 'public', SetAccess = 'protected')
        m_fx; % functional value, a scalar value or a vector of scalars corresponding to the number of scenarios
        m_fx0; % initial functional value
        m_dfdx; % gradient of functional, a matrix of size (ny x nx) or (ny x nx x numScenarios)
        m_solver; % solver
        m_numEvaluations; % number of evaluations
        m_isConstraint; % true if the functional is constraint, false if the functional is objective
        m_upperBound; % upper bound
        m_scale;% scaling factor for gradients
    end
    methods (Abstract)
        % abstract methods
        % evaluate the functional value at design and state variables queried from the solver
        % input: obj, ubound (optional)
        % output: obj, functional value
        % ubound is optional, if not provided, the functional is considered as objective and returns fx
        % if ubound is provided, the functional is considered as constraint and returns gx := fx/ubound - 1 <= 0
        [obj,value] = evaluate(obj)

        % gradient of the functional at design and state variables queried from the solver.
        % additional adjoint problems are solved in this function to obtain the lagrnage multipliers
        % input: obj, ubound (optional)
        % output: obj, functional gradient w.r.t. design variables
        % ubound is optional, if not provided, the functional is considered as objective and returns df/dx
        % if ubound is provided, the functional is considered as constraint and returns dgx/dx := dfx/dx/ubound
        [obj,grad] = gradient(obj)
    end
    methods
        %% CONSTRUCTOR
        function obj = functional(solver, ubound)
            % constructor
            if (~isa(solver, 'simulation2d')), error('solver must be an instance of simulation2d class!');end % check if solver is valid

            obj.m_isConstraint = false; % objective
            if (~isnan(ubound))
                obj.m_isConstraint = true; % constraint
                obj.m_upperBound = ubound;
            end
            obj.m_solver = solver; % assign solver
            obj.m_numEvaluations = 0; % initialize number of evaluations
            obj.m_scale = 1.0;
        end
    end
end
