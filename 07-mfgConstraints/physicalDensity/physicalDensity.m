%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discription:                                                              %
% This is a class for imposing physical density constraints used in     %
% shape and topology optimization (e.g., SIMP).                             %
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

classdef  physicalDensity < mfgConstraints
    properties (GetAccess = 'public', SetAccess = 'protected')
        m_beta; % Sharpness of Heaviside projection
        m_eta; % Threshold of Heaviside projection

        m_betaMax;
        % number of iterations for optimization
        % updates based on sensitivity filtering
        m_numIter;
    end
    methods
        %% CONSTRUCTOR
        function obj = physicalDensity(solver,beta,eta,betaMax)
            % check if solver is valid
            if (~isa(solver, 'simulation2d')), error('solver must be an instance of simulation2d class!');end % check if solver is valid
            
            if (nargin < 2)
                beta = 1;
            end
            if (nargin < 3)
                eta = 0.5;
            end
             if (nargin < 4)
                betaMax = 4;
            end
            % constructor based on superclass
            obj = obj@mfgConstraints(solver);
            obj.m_beta = beta;
            obj.m_eta = eta;
            obj.m_numIter = 0;
            obj.m_betaMax = betaMax;
        end

        function obj = setParameters(obj,beta,eta,betaMax)
            % Reset the parameters for a new optimization run
            obj.m_beta = beta;
            obj.m_eta = eta;
            obj.m_betaMax = betaMax;
        end 

        function [filteredDesign] = filterDesign(obj, design)
            % filter the design variables
            % input: obj, design variables
            % output: obj, filtered design variables
            numerator = tanh(obj.m_beta*obj.m_eta) + tanh(obj.m_beta*(design - obj.m_eta));
            denominator = tanh(obj.m_beta*obj.m_eta) + tanh(obj.m_beta*(1 - obj.m_eta));
            filteredDesign = numerator/denominator;
            filteredDesign = filteredDesign .* obj.m_solver.m_existingElems;
        end

        function [filteredSensitivity] = filterSensitivity(obj, design, sensField)
            % filter the sensitivity fields
            % input: obj, sensitivity fields
            % output: obj, filtered sensitivity fields
            denominator = tanh(obj.m_beta*obj.m_eta) + tanh(obj.m_beta*(1 - obj.m_eta));

            gradFilteredDesign = obj.m_beta * (1 - tanh(obj.m_beta*(design - obj.m_eta)).^2) / denominator;

            filteredSensitivity = gradFilteredDesign .* sensField;
            filteredSensitivity = filteredSensitivity .* obj.m_solver.m_existingElems;

            obj.m_numIter = obj.m_numIter + 1;

            % update beta every 50 iterations to increase sharpness
            if (mod(obj.m_numIter, 50) == 0)
                obj.m_beta = min(obj.m_beta * 2, obj.m_betaMax);
            end
        end

    end
end
