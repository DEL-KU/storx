%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                              %
% Random Search parameter optimization (concrete subclass).                 %
% Samples random feasible points in parameter space and tracks the best.   %
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

classdef parameterOpt2d_RS < parameterOpt2d

    properties (GetAccess = 'public', SetAccess = 'protected')
        % Number of feasible random samples to collect before stopping
        m_numRandomSamples = 10;

        % Maximum number of total objective evaluations
        m_numRandomFuncCounts = 100;
    end

    methods

        function obj = parameterOpt2d_RS(brepHandle, solverHandle, param, ...
                objective, constraints, exportGIF, testMode)

            if nargin < 6 || isempty(exportGIF), exportGIF = false; end
            if nargin < 7 || isempty(testMode),  testMode  = false; end

            obj = obj@parameterOpt2d(brepHandle, solverHandle, param, ...
                objective, constraints, [], [], exportGIF, testMode);
        end

        function obj = setNumberOfRandomSearchSamples(obj, nSamples, nFuncCounts)
            obj.m_numRandomSamples = nSamples;
            if nargin >= 3
                obj.m_numRandomFuncCounts = nFuncCounts;
            end
        end

        function obj = optimize(obj)
            obj.openProgressBar('RS: Random Search');
            LB = obj.m_param0.lb ./ obj.m_param0.value;
            UB = obj.m_param0.ub ./ obj.m_param0.value;
            x0 = ones(size(LB));

            nSamples   = obj.m_numRandomSamples;
            nFuncCount = obj.m_numRandomFuncCounts;
            obj.m_totalEvals = nFuncCount;
            bestObj    = Inf;
            bestX      = [];

            solutionsRS  = [];
            funcCount    = 0;
            feasSolCount = 0;

            while feasSolCount < nSamples && funcCount < nFuncCount
                xCandidate = LB + rand(size(LB)) .* (UB - LB);

                [c, ceq] = obj.evaluateConstraint(xCandidate);
                fval     = obj.evaluateObjective(xCandidate);
                funcCount = funcCount + 1;

                isFeasible = all(c <= 0) && all(abs(ceq) <= 1e-6);
                if isFeasible
                    feasSolCount = feasSolCount + 1;
                    solutionsRS(feasSolCount).X    = xCandidate; %#ok
                    solutionsRS(feasSolCount).Fval = fval;       %#ok
                    if fval < bestObj
                        bestObj = fval;
                        bestX   = xCandidate;
                    end
                end
            end

            if ~isempty(bestX)
                xMin = bestX;
                flag = 1;
            else
                xMin = x0;
                flag = -1;
            end

            obj.m_feasibleExploredSolutions = solutionsRS;
            obj = obj.finalizeResults(xMin, flag, funcCount);

            if ~obj.m_testMode
                obj = obj.plotSampleSpace();
            end
        end

    end
end
