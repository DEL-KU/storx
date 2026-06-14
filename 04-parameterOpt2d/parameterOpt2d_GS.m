%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                              %
% Global Search parameter optimization (concrete subclass).                 %
% Uses MATLAB's GlobalSearch with fmincon to escape local minima.          %
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

classdef parameterOpt2d_GS < parameterOpt2d

    methods

        function obj = parameterOpt2d_GS(brepHandle, solverHandle, param, ...
                objective, constraints, ...
                terminationTolerance, finiteDifferenceStepSize, ...
                exportGIF, testMode)

            if nargin < 8 || isempty(exportGIF), exportGIF = false; end
            if nargin < 9 || isempty(testMode),  testMode  = false; end

            obj = obj@parameterOpt2d(brepHandle, solverHandle, param, ...
                objective, constraints, ...
                terminationTolerance, finiteDifferenceStepSize, ...
                exportGIF, testMode);
        end

        function obj = optimize(obj)
            obj.openProgressBar('GS: Global Search');
            LB = obj.m_param0.lb ./ obj.m_param0.value;
            UB = obj.m_param0.ub ./ obj.m_param0.value;
            x0 = ones(size(LB));

            opt = optimoptions('fmincon', ...
                'TolX',                     obj.m_terminationTolerance, ...
                'TolFun',                   obj.m_terminationTolerance, ...
                'ConstraintTolerance',      obj.m_terminationTolerance, ...
                'FiniteDifferenceStepSize', obj.m_finiteDifferenceStepSize);

            problem = createOptimProblem('fmincon', ...
                'objective', @obj.evaluateObjective, ...
                'x0', x0, 'lb', LB, 'ub', UB, ...
                'nonlcon', @obj.evaluateConstraint, 'options', opt);

            gs = GlobalSearch('MaxTime', 300, 'NumTrialPoints', 300, 'NumStageOnePoints', 100);
            obj.m_totalEvals = gs.NumTrialPoints * (1 + obj.m_numParams);
            [xMin, ~, flag, output, solutions] = run(gs, problem);

            obj.m_feasibleExploredSolutions = solutions;
            funcCount = output.funcCount * (1 + obj.m_numParams);
            obj = obj.finalizeResults(xMin, flag, funcCount);

            if ~obj.m_testMode
                obj = obj.plotSampleSpace();
            end
        end

    end
end
