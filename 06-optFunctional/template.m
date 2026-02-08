%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % %           You can use this as a template for your own class.      % % %
% Discription:                                                              %
% This is a template class for functionals used in shape and                %
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

classdef template < functional %% change this to the name of your class
    properties (GetAccess = 'public', SetAccess = 'protected')        
        %% ADD YOUR OWN PROPERTIES SPECIFIC TO YOUR FUNCTIONAL
        %% ...
    end

    methods
        %% CONSTRUCTOR
        function obj = template(solver, ub) % change this to the name of your class (should match the file name)
            %   input: solver, ub (optional)
            %   output: obj
            %   ub is optional, if not provided, the functional is considered as objective and returns fx
            %   if ub is provided, the functional is considered as constraint and returns gx := fx/ubound - 1 <= 0

            % check if solver is valid (change SOLVER_NAME to the name of your class)
            %  for example:
            % if (~isa(solver, 'fea2d_elasticity')), error('solver must be an instance of fea2d_elasticity class!');end

            if (~isa(solver, 'SOLVER_NAME')), error('solver must be an instance of SOLVER_NAME class!');end
           
            % if upper bound is provided set the value, otherwise set it to NaN
            upper_bound = NaN;
            if (nargin > 1)
                upper_bound = ub;
            end

            % constructor based on superclass functional
            obj = obj@functional(solver, upper_bound);
        end

        % evaluate the YOUR_FUNCTIONAL value at design and state variables queried from the SOLVER_NAME solver
        % input: obj
        % output: obj, YOUR_FUNCTIONAL_VALUE
        % if objective, the functional is considered as objective and returns fx
        % otherwise, the functional is considered as constraint and returns gx := fx/ubound - 1 <= 0
        function [obj,value] = evaluate(obj)
            %% ADD YOUR OWN CODE HERE


            %% ...


            %% DO NOT CHANGE BEYON THIS LINE
            if (obj.m_numEvaluations == 0), obj.m_fx0 = obj.m_fx; end
                
            if (obj.m_isConstraint)
                value = obj.m_fx/obj.m_upperBound - 1; % as constraint
            else
                value = obj.m_fx; % as objective
            end
            obj.m_numEvaluations = obj.m_numEvaluations + 1; % update number of evaluations
        end

        % gradient of the YOUR_FUNCTIONAL at design and state variables queried from the SOLVER_NAME solver.
        % additional adjoint problems are solved in this function to obtain the lagrnage multipliers
        % input: obj, ubound (optional)
        % output: obj, YOUR_FUNCTIONAL_GRADIENT w.r.t. design variables
        % if objective, the functional is considered as objective and returns df/dx
        % otherwise, the functional is considered as constraint and returns dgx/dx := dfx/dx/ubound
        function [obj,grad] = gradient(obj)
            %% ADD YOUR OWN CODE HERE


            %% ...

            %% DO NOT CHANGE BEYON THIS LINE
            if (obj.m_isConstraint)
                grad = obj.m_dfdx/obj.m_upperBound;
            else
                grad = obj.m_dfdx;
            end

        end
    end
    
    methods (Access = 'private')
        %% ADD YOUR OWN PRIVATE METHODS
        % These methods are not called from outside the class
        % and can be used to implement special functions like solving adjoint problems
        % or other helper functions to keep evaluate and gradient functions clean
        function obj = some_private_method(obj)
             %% ADD YOUR OWN CODE HERE


            %% ...
        end

    end
end
