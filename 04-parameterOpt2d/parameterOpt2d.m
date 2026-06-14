%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                              %
% Abstract base class for 2D parameter-based geometry optimization.        %
% Defines shared infrastructure (solver management, objective/constraint   %
% evaluation, history tracking, and visualization) used by all concrete    %
% method subclasses: parameterOpt2d_FD, parameterOpt2d_RS,                 %
% parameterOpt2d_GS, parameterOpt2d_MS.                                    %
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

classdef (Abstract) parameterOpt2d < handle

    properties (GetAccess = 'public', SetAccess = 'protected')
        % Function handle to create the geometry (brep) from parameters
        m_brepHandle

        % Function handle to build and return the FEA solver object
        m_solverHandle

        % The initial and final solver objects (contain geometry, FEA data, etc.)
        m_solverInitial
        m_solverFinal

        % Parameter specification struct with fields:
        %   .value : nominal parameter values
        %   .lb    : lower bound
        %   .ub    : upper bound
        m_param0

        % Computed optimal parameter values after optimization
        m_optimalParams

        % Number of design parameters
        m_numParams

        % Objective type (e.g. 'compliance')
        m_objective

        % Constraints struct with possible fields:
        %   .area      : area constraint value
        %   .perimeter : perimeter constraint value
        %   .type      : 'ineq' or 'eq'
        m_constraints

        % Target area or perimeter extracted from constraints
        m_targetArea
        m_targetPerimeter

        % Baseline performance metrics for the initial design
        m_cx0
        m_area0
        m_perim0
        m_maxDef0
        m_maxStress0

        % Termination tolerance for fmincon
        m_terminationTolerance

        % Finite difference step size for fmincon
        m_finiteDifferenceStepSize

        % Flag for saving .gif animations during optimization
        m_exportGIF

        % Structure of solutions found by certain methods (e.g. GlobalSearch)
        m_feasibleExploredSolutions

        % Unified results structure: stores both final results and history
        m_results

        % Flag for test mode - will not plot or export GIFs
        m_testMode = false;
    end

    properties (Access = protected)
        % Iteration history stored on the object (replaces persistent variables)
        m_paramsHistory
        m_objectiveHistory
        m_feasibleHistory
        m_areaHistory
        m_perimHistory

        % Waitbar handle and total eval count for progress display
        m_progressBar  = []
        m_totalEvals   = Inf;
    end

    methods (Abstract)
        obj = optimize(obj)
    end

    methods

        %% Constructor
        function obj = parameterOpt2d(brepHandle, solverHandle, param, ...
                objective, constraints, ...
                terminationTolerance, finiteDifferenceStepSize, ...
                exportGIF, testMode)
            % Constructor for parameterOpt2d base class.
            %  Inputs:
            %   brepHandle             : function handle to create geometry from parameters
            %   solverHandle           : function handle to create FEA solver object
            %   param                  : struct with .value .lb .ub
            %   objective              : string, e.g. 'compliance'
            %   constraints            : struct with fields like .area, .perimeter, .type
            %   terminationTolerance   : termination tol for fmincon (optional)
            %   finiteDifferenceStepSize : step size for finite diff in fmincon (optional)
            %   exportGIF              : boolean to export .gif files (optional)
            %   testMode               : boolean, suppresses plotting/export (optional)

            if nargin < 6 || isempty(terminationTolerance)
                terminationTolerance = 1e-6;
            end
            if nargin < 7 || isempty(finiteDifferenceStepSize)
                finiteDifferenceStepSize = 1e-6;
            end
            if nargin < 8 || isempty(exportGIF)
                exportGIF = false;
            end
            if nargin < 9 || isempty(testMode)
                testMode = false;
            end

            obj.m_brepHandle              = brepHandle;
            obj.m_solverHandle            = solverHandle;
            obj.m_param0                  = param;
            obj.m_numParams               = numel(param.value);
            obj.m_objective               = objective;
            obj.m_constraints             = constraints;
            obj.m_terminationTolerance    = terminationTolerance;
            obj.m_finiteDifferenceStepSize = finiteDifferenceStepSize;
            obj.m_exportGIF               = exportGIF;
            obj.m_testMode                = testMode;
            if testMode, obj.m_exportGIF = false; end

            % Validate solver type
            brep0 = obj.m_brepHandle(param.value);
            obj.m_solverInitial = obj.m_solverHandle(brep0);
            if ~isa(obj.m_solverInitial, 'fea2d_elasticity') && ...
                    ~isa(obj.m_solverInitial, 'triFEA2d_elasticity')
                error('parameterOpt2d: solver must be fea2d_elasticity or triFEA2d_elasticity.');
            end

            % Solve the initial geometry to get baseline metrics
            scenarioId = 1;
            obj.m_solverInitial  = obj.solve(param.value);
            obj.m_cx0            = obj.m_solverInitial.computeCompliance();
            obj.m_maxDef0        = obj.m_solverInitial.m_maxDef(scenarioId);
            obj.m_maxStress0     = obj.m_solverInitial.m_maxStress(scenarioId);
            obj.m_area0          = obj.m_solverInitial.getArea();
            obj.m_perim0         = obj.m_solverInitial.getPerimeter();

            % Initialize history arrays
            obj.m_paramsHistory    = [];
            obj.m_objectiveHistory = [];
            obj.m_feasibleHistory  = [];
            obj.m_areaHistory      = [];
            obj.m_perimHistory     = [];

            % Initialize unified results struct
            obj.m_results = struct( ...
                'flag',              [], ...
                'initialParams',     [], ...
                'initialArea',       [], ...
                'initialPerimeter',  [], ...
                'initialCompliance', [], ...
                'initialDef',        [], ...
                'initialStress',     [], ...
                'finalParams',       [], ...
                'finalArea',         [], ...
                'finalPerimeter',    [], ...
                'finalCompliance',   [], ...
                'finalDef',          [], ...
                'finalStress',       [], ...
                'nFEAs',             [], ...
                'history', struct( ...
                    'x',        [], ...
                    'f',        [], ...
                    'feasible', [], ...
                    'area',     [], ...
                    'perim',    [] ) );

            % Extract constraint targets
            if isfield(obj.m_constraints, 'area')
                obj.m_targetArea = obj.m_constraints.area;
            elseif isfield(obj.m_constraints, 'perimeter')
                obj.m_targetPerimeter = obj.m_constraints.perimeter;
            else
                disp('No constraint imposed!')
            end

            if ~obj.m_testMode
                obj.m_solverInitial.plotBoundaryCondition();
                obj.m_solverInitial.plotGeometry(1, 0, 'Initial Geometry');
            end
        end

        %% Solver creation and execution
        function fem = solve(obj, params)
            fem = obj.createSolver(params);
            fem = fem.preProcess();
            fem = fem.solve();
            fem = fem.postProcess();

            if ~obj.m_testMode
                plt = PlotId;
                fem.plotGeometry(plt.brep_optimized, 0, 'Optimized Geometry');
                fem.printElascticityResults();
            end
        end

        function fem = createSolver(obj, params)
            brep = obj.m_brepHandle(params);
            fem  = obj.m_solverHandle(brep);
        end

        function fem = solveQuiet(obj, params)
            fem = obj.createSolver(params);
            fem = fem.preProcess();
            fem = fem.solve();
            fem = fem.postProcess();
        end

        %% OBJECTIVE & CONSTRAINTS
        function fx = evaluateObjective(obj, x)
            params = x .* obj.m_param0.value;

            if strcmp(obj.m_objective, 'compliance')
                fem = obj.solve(params);
                cx  = fem.computeCompliance();
                fx  = cx / obj.m_cx0;
            else
                disp(['Objective ' obj.m_objective ' is not implemented!']);
                fx = 0;
                cx = 0;
            end

            obj.m_paramsHistory    = [obj.m_paramsHistory;    params];
            obj.m_objectiveHistory = [obj.m_objectiveHistory; cx];
            inBounds = all(params <= obj.m_param0.ub & params >= obj.m_param0.lb);
            obj.m_feasibleHistory  = [obj.m_feasibleHistory;  double(inBounds)];

            if ~obj.m_testMode
                n = numel(obj.m_objectiveHistory);
                fprintf('  eval %-4d | C = %.6f\n', n, cx);
            end
            obj.updateProgressBar(cx);
            if obj.m_exportGIF, export_gifs(); end
        end

        function [cineq, ceq] = evaluateConstraint(obj, x)
            if isfield(obj.m_constraints, 'area')
                [cineq, ceq] = obj.evaluateAreaConstraint(x);
            elseif isfield(obj.m_constraints, 'perimeter')
                [cineq, ceq] = obj.evaluatePerimeterConstraint(x);
            else
                disp('No constraint imposed!')
                cineq = [];
                ceq   = [];
            end
        end

        function [cineq, ceq] = evaluateAreaConstraint(obj, x)
            params = x .* obj.m_param0.value;
            fem    = obj.createSolver(params);
            area   = fem.getArea();

            if strcmp(obj.m_constraints.type, 'ineq')
                cineq = area / obj.m_targetArea - 1.0;
                ceq   = [];
            else
                cineq = [];
                ceq   = area / obj.m_targetArea - 1.0;
            end

            obj.m_areaHistory = [obj.m_areaHistory; area];
        end

        function [cineq, ceq] = evaluatePerimeterConstraint(obj, x)
            params = x .* obj.m_param0.value;
            fem    = obj.createSolver(params);
            perim  = fem.getPerimeter();

            if strcmp(obj.m_constraints.type, 'ineq')
                cineq = perim / obj.m_targetPerimeter - 1.0;
                ceq   = [];
            else
                cineq = [];
                ceq   = perim / obj.m_targetPerimeter - 1.0;
            end

            obj.m_perimHistory = [obj.m_perimHistory; perim];
        end

        %% POST-OPTIMIZATION RESULTS
        function obj = finalizeResults(obj, xMin, flag, funcCount)
            % finalizeResults
            % Called by each concrete optimize() to store the shared post-opt results.

            scenarioId = 1;
            obj.m_results.flag = flag;

            obj.m_results.initialParams     = obj.m_param0.value;
            obj.m_results.initialArea       = obj.m_area0;
            obj.m_results.initialPerimeter  = obj.m_perim0;
            obj.m_results.initialCompliance = obj.m_cx0;
            obj.m_results.initialDef        = obj.m_maxDef0;
            obj.m_results.initialStress     = obj.m_maxStress0;

            obj.m_optimalParams           = xMin .* obj.m_param0.value;
            fem = obj.solve(obj.m_optimalParams);

            obj.m_results.finalParams     = obj.m_optimalParams;
            obj.m_results.finalArea       = fem.getArea();
            obj.m_results.finalPerimeter  = fem.getPerimeter();
            obj.m_results.finalCompliance = fem.computeCompliance();
            obj.m_results.finalDef        = fem.m_maxDef(scenarioId);
            obj.m_results.finalStress     = fem.m_maxStress(scenarioId);
            obj.m_results.nFEAs           = funcCount;
            obj.m_solverFinal             = fem;

            obj.m_results.history.x        = obj.m_paramsHistory;
            obj.m_results.history.f        = obj.m_objectiveHistory;
            obj.m_results.history.feasible = obj.m_feasibleHistory;
            obj.m_results.history.area     = obj.m_areaHistory;
            obj.m_results.history.perim    = obj.m_perimHistory;

            obj.closeProgressBar();
            if ~obj.m_testMode
                disp('Exit Flag: ');         disp(obj.m_results.flag);
                disp('Initial Params: ');    disp(obj.m_results.initialParams);
                disp('Final Params: ');      disp(obj.m_results.finalParams);
                disp('Initial/Final Area: ');
                disp([obj.m_results.initialArea       obj.m_results.finalArea]);
                disp('Initial/Final Perimeter: ');
                disp([obj.m_results.initialPerimeter  obj.m_results.finalPerimeter]);
                disp('Initial/Final Compliance: ');
                disp([obj.m_results.initialCompliance obj.m_results.finalCompliance]);
                disp('Initial/Final Def: ');
                disp([obj.m_results.initialDef        obj.m_results.finalDef]);
                disp('Initial/Final Stress: ');
                disp([obj.m_results.initialStress     obj.m_results.finalStress]);
                disp('#FEAs:'); disp(obj.m_results.nFEAs);
            end
        end

        %% PROGRESS BAR
        function openProgressBar(obj, title)
            if obj.m_testMode, return; end
            obj.m_progressBar = waitbar(0, 'Starting...', 'Name', title);
        end

        function updateProgressBar(obj, cx)
            if obj.m_testMode || isempty(obj.m_progressBar) || ~isvalid(obj.m_progressBar)
                return;
            end
            n = numel(obj.m_objectiveHistory);
            if isfinite(obj.m_totalEvals) && obj.m_totalEvals > 0
                frac = min(n / obj.m_totalEvals, 1);
                msg  = sprintf('Eval %d / %d  |  C = %.4f', n, obj.m_totalEvals, cx);
            else
                frac = min((n - 1) / 150, 0.99);
                msg  = sprintf('Eval #%d  |  C = %.4f', n, cx);
            end
            waitbar(frac, obj.m_progressBar, msg);
        end

        function closeProgressBar(obj)
            if ~isempty(obj.m_progressBar) && isvalid(obj.m_progressBar)
                close(obj.m_progressBar);
                obj.m_progressBar = [];
            end
        end

        %% VISUALIZATION
        function obj = plotConvergence(obj)
            plt = PlotId;
            figure(plt.convergence);
            set(gcf, 'Name', 'Convergence');

            nIter = numel(obj.m_results.history.f);

            if isfield(obj.m_constraints, 'area')
                consVals  = obj.m_results.history.area(1:nIter);
                consLabel = 'Area ($A$)';
            elseif isfield(obj.m_constraints, 'perimeter')
                consVals  = obj.m_results.history.perim(1:nIter);
                consLabel = 'Perimeter ($P$)';
            else
                disp('No constraint imposed!');
                consVals  = zeros(1, nIter);
                consLabel = 'No Constraint';
            end

            f = obj.m_results.history.f(1:nIter);
            plot(1:nIter, f, '-b', 'LineWidth', 2);
            hold on; grid on;
            xlabel('Iteration');

            if strcmp(obj.m_objective, 'compliance')
                ylabel('Compliance $C$');
            else
                ylabel('Objective $\varphi$');
            end
            axis tight;

            yyaxis right
            plot(1:nIter, consVals, '--r', 'LineWidth', 2);
            ylabel(consLabel);
            hold off;
        end

        function obj = plotSampleSpace(obj)
            plt = PlotId;

            useArea      = isfield(obj.m_constraints, 'area');
            usePerimeter = isfield(obj.m_constraints, 'perimeter');

            params = obj.m_results.history.x;
            fVals  = obj.m_results.history.f(:);
            nSamples = min(size(params, 1), numel(fVals));
            if nSamples == 0
                warning('No samples recorded in optimization history!');
                return;
            end

            params = params(1:nSamples, :);
            fVals  = fVals(1:nSamples);
            sampleIds = (1:nSamples)';
            xVals = sampleIds;
            xLabel = 'Sample number';
            targetLabel = '';
            targetValue = [];
            hasConstraintAxis = useArea || usePerimeter;

            if hasConstraintAxis
                xVals = zeros(nSamples, 1);
                if useArea
                    xLabel = 'Area, $A$';
                    targetLabel = '$A_{\mathrm{target}}$';
                    targetValue = obj.m_targetArea;
                elseif usePerimeter
                    xLabel = 'Perimeter, $P$';
                    targetLabel = '$P_{\mathrm{target}}$';
                    targetValue = obj.m_targetPerimeter;
                end

                for i = 1:nSamples
                    xVals(i) = obj.getParamMetricForPlot(params(i, :), useArea);
                end
            end

            finiteSamples = isfinite(xVals) & isfinite(fVals);
            if hasConstraintAxis
                tol = max(obj.m_terminationTolerance, 1e-6) * max(1, abs(targetValue));
                if isfield(obj.m_constraints, 'type') && strcmp(obj.m_constraints.type, 'eq')
                    feasibleSamples = abs(xVals - targetValue) <= tol;
                else
                    feasibleSamples = xVals <= targetValue + tol;
                end
            else
                feasibleSamples = true(nSamples, 1);
            end
            feasibleSamples = feasibleSamples & finiteSamples;
            infeasibleSamples = finiteSamples & ~feasibleSamples;

            figure(plt.pareto_front);
            clf;
            set(gcf, 'Name', 'SampleSpace');
            ax = gca;
            ax.TickLabelInterpreter = 'latex';
            hold on; grid on;

            if any(infeasibleSamples)
                scatter(xVals(infeasibleSamples), fVals(infeasibleSamples), ...
                    32, [0.70 0.70 0.70], 'o', 'DisplayName', 'Infeasible samples');
            end
            if any(feasibleSamples)
                scatter(xVals(feasibleSamples), fVals(feasibleSamples), ...
                    38, sampleIds(feasibleSamples), 'filled', ...
                    'MarkerEdgeColor', [0.15 0.15 0.15], ...
                    'DisplayName', 'Feasible samples');
                cb = colorbar;
                cb.TickLabelInterpreter = 'latex';
                cb.Label.Interpreter = 'latex';
                cb.Label.String = 'Evaluation number';
            end

            if nnz(feasibleSamples) >= 2
                feasibleIds = find(feasibleSamples);
                [~, order] = sort(xVals(feasibleIds));
                envelopeIds = [];
                bestSoFar = Inf;
                for k = 1:numel(order)
                    id = feasibleIds(order(k));
                    if fVals(id) < bestSoFar
                        envelopeIds(end+1) = id; %#ok<AGROW>
                        bestSoFar = fVals(id);
                    end
                end
                if numel(envelopeIds) >= 2
                    plot(xVals(envelopeIds), fVals(envelopeIds), '-k', ...
                        'LineWidth', 1.5, 'DisplayName', 'Sample envelope');
                end
            end

            [xInitial, xFinal] = obj.getDesignMetricForPlot(hasConstraintAxis, useArea);
            scatter(xInitial, obj.m_results.initialCompliance, 90, 's', ...
                'MarkerFaceColor', [0.00 0.45 0.74], ...
                'MarkerEdgeColor', 'k', 'DisplayName', 'Initial design');
            scatter(xFinal, obj.m_results.finalCompliance, 110, 'p', ...
                'MarkerFaceColor', [0.85 0.10 0.10], ...
                'MarkerEdgeColor', 'k', 'DisplayName', 'Final design');

            [xSolution, ySolution] = obj.getExploredSolutionMetricsForPlot(hasConstraintAxis, useArea);
            if ~isempty(xSolution)
                scatter(xSolution, ySolution, 70, 'd', ...
                    'MarkerEdgeColor', [0.49 0.18 0.56], ...
                    'LineWidth', 1.2, 'DisplayName', 'Local solution endpoints');
            end

            if hasConstraintAxis
                yLimits = ylim;
                plot([targetValue targetValue], yLimits, '--', ...
                    'Color', [0.20 0.20 0.20], 'LineWidth', 1.2, ...
                    'DisplayName', targetLabel);
                ylim(yLimits);
            end

            xlabel(xLabel, 'Interpreter', 'latex');
            if strcmp(obj.m_objective, 'compliance')
                ylabel('Compliance, $C$', 'Interpreter', 'latex');
            else
                ylabel('Objective, $\varphi$', 'Interpreter', 'latex');
            end
            legend('Location', 'best', 'Interpreter', 'latex');
            axis tight;
            hold off;

            fprintf('Sample space: %d evaluations, %d feasible\n', ...
                nnz(finiteSamples), nnz(feasibleSamples));
        end

    end

    methods (Access = protected)

        function [xInitial, xFinal] = getDesignMetricForPlot(obj, hasConstraintAxis, useArea)
            if hasConstraintAxis
                if useArea
                    xInitial = obj.m_results.initialArea;
                    xFinal   = obj.m_results.finalArea;
                else
                    xInitial = obj.m_results.initialPerimeter;
                    xFinal   = obj.m_results.finalPerimeter;
                end
            else
                xInitial = 1;
                xFinal   = numel(obj.m_results.history.f);
            end
        end

        function [xSolution, ySolution] = getExploredSolutionMetricsForPlot(obj, hasConstraintAxis, useArea)
            xSolution = [];
            ySolution = [];
            if ~hasConstraintAxis || isempty(obj.m_feasibleExploredSolutions)
                return;
            end

            solutions = obj.m_feasibleExploredSolutions;
            for i = 1:numel(solutions)
                if isempty(solutions(i).X) || isempty(solutions(i).Fval)
                    continue;
                end

                paramCandidate = solutions(i).X .* obj.m_param0.value;
                xCandidate = obj.getParamMetricForPlot(paramCandidate, useArea);
                yCandidate = solutions(i).Fval * obj.m_results.initialCompliance;

                xSolution(end+1, 1) = xCandidate; %#ok<AGROW>
                ySolution(end+1, 1) = yCandidate; %#ok<AGROW>
            end
        end

        function metricValue = getParamMetricForPlot(obj, params, useArea)
            brepCandidate = brep2d(obj.m_brepHandle(params));
            if useArea
                metricValue = brepCandidate.getArea();
            else
                metricValue = brepCandidate.getPerimeter();
            end
        end

    end
end
