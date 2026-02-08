%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                              %
% This class implements a 2D parameter-based geometry optimization using    %
% different search methods (FD, RS, GS, MS). The primary goal is to         %
% modify design parameters to minimize an objective function (e.g.          %
% compliance) under constraints (e.g. on area or perimeter).                %
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

classdef parameterOpt2d
    % parameterOpt2d
    % This class implements a 2D parameter-based geometry optimization
    % using different search methods (FD, RS, GS, MS). The primary goal
    % is to modify design parameters to minimize an objective function
    % (e.g. compliance) under constraints (e.g. on area or perimeter).

    properties (GetAccess = 'public', SetAccess = 'protected')
        % Optimization method:
        %  - RS: Random Search
        %  - FD: Finite Difference
        %  - MS: Multi-Start
        %  - GS: Global Search
        m_method (1,1) string = "RS"

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

        % Number of random feasible samples explored
        m_numRandomSamples;

        % Number of random samples evaluated
        m_numRandomFuncCounts;

        % Number of local fmincon problems from each multi-start starting
        % point
        m_numLocalMultiStartProblems;

        % Flag for saving .gif animations during optimization
        m_exportGIF

        % Structure of solutions found by certain methods (e.g. GlobalSearch)
        m_feasibleExploredSolutions

        % Unified results structure: stores both final results and history
        m_results

        % Flag for test mode, will not plot or export GIFs
        m_testMode = false;
    end

    methods

        %% Constructor
        function obj = parameterOpt2d(brepHandle, solverHandle, param, ...
                objective, constraints, ...
                terminationTolerance, finiteDifferenceStepSize, ...
                method, exportGIF, testMode)
            % Constructor for parameterOpt2d class
            %  Inputs:
            %   brepHandle  : function handle to create geometry from parameters
            %   solverHandle: function handle to create FEA solver object
            %   param       : struct with .value .lb .ub
            %   objective   : string representing objective function, e.g. 'compliance'
            %   constraints : struct with fields like .area, .perimeter, .type
            %   terminationTolerance: termination tol for fmincon
            %   finiteDifferenceStepSize: step size for finite diff in fmincon
            %   method      : search method (FD, RS, GS, MS)
            %   exportGIF     : boolean to export .gif files


            % default values for optional parameters
            if nargin < 6 || isempty(terminationTolerance)
                terminationTolerance = 1e-6; % default tolerance
            end
            if nargin < 7 || isempty(finiteDifferenceStepSize)
                finiteDifferenceStepSize = 1e-6; % default step size
            end
            if nargin < 8 || isempty(method)
                method = "RS"; % default method
            end
            if nargin < 9 || isempty(exportGIF)
                exportGIF = false;
            end
            if nargin < 10 || isempty(testMode)
                testMode = false;
            end

            % Store function handles and input configuration
            obj.m_brepHandle = brepHandle;
            obj.m_solverHandle = solverHandle;
            obj.m_param0 = param;
            obj.m_numParams = numel(param.value);
            obj.m_objective = objective;
            obj.m_constraints = constraints;
            obj.m_terminationTolerance = terminationTolerance;
            obj.m_finiteDifferenceStepSize = finiteDifferenceStepSize;
            obj.m_method = method;
            obj.m_exportGIF = exportGIF;
            obj.m_testMode = testMode;
            if testMode, obj.m_exportGIF = false; end
            % Create initial design and compute its performance
            brep0 = obj.m_brepHandle(param.value);
            obj.m_solverInitial = obj.m_solverHandle(brep0);

            % check the solverHandle
            if ~isa(obj.m_solverInitial, 'fea2d_elasticity') && ~isa(obj.m_solverInitial, 'triFEA2d_elasticity')
                error('solver must be an instance of fea2d_elasticity class!');
            end

            % Solve the initial geometry to get baseline results
            scenarioId = 1;
            obj.m_solverInitial = obj.solve(param.value);
            obj.m_cx0       = obj.m_solverInitial.computeCompliance();
            obj.m_maxDef0   = obj.m_solverInitial.m_maxDef(scenarioId);
            obj.m_maxStress0= obj.m_solverInitial.m_maxStress(scenarioId);
            obj.m_area0     = obj.m_solverInitial.getArea();
            obj.m_perim0    = obj.m_solverInitial.getPerimeter();

            % Initialize the unified results struct
            obj.m_results = struct( ...
                'flag',      [], ...
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
                'history',           struct( ...
                'x',          [], ...
                'f',          [], ...
                'feasible',   [], ...
                'area',       [], ...
                'perim',      [] ) ...
                );

            % If an area or perimeter constraint is provided, store it
            if isfield(obj.m_constraints, 'area')
                obj.m_targetArea = obj.m_constraints.area;
            elseif isfield(obj.m_constraints, 'perimeter')
                obj.m_targetPerimeter = obj.m_constraints.perimeter;
            else
                disp('No constraint imposed!')
            end

            % Prime persistent variables in objective/constraint methods
            obj.evaluateObjective([],1);
            obj.evaluateConstraint([],1);

            if ~obj.m_testMode
                % Plot initial boundary conditions and geometry
                obj.m_solverInitial.plotBoundaryCondition();
                obj.m_solverInitial.plotGeometry(1,0,'Initial Geometry');
            end
        end

        %% Set number of random samples
        function obj = setNumberOfRandomSearchSamples(obj,nSamples,nFuncCounts)
            obj.m_numRandomSamples = nSamples;
            obj.m_numRandomFuncCounts = nFuncCounts;
        end

        %% Set number of multi-start local problems
        function obj = setNumberOfMultiStartLocalProblems(obj,nLocalProblems)
            obj.m_numLocalMultiStartProblems = nLocalProblems;
        end

        %% Solver creation and execution
        function fem = solve(obj, params)
            % solve
            % Creates a solver with given params, runs FEA, and post-processes
            fem = obj.createSolver(params);
            fem = fem.preProcess();  % FEA pre-processing
            fem = fem.solve();       % Solve FEA
            fem = fem.postProcess(); % Post-processing of results

            if ~obj.m_testMode
                % Optionally plot the updated geometry
                plt = PlotId;
                fem.plotGeometry(plt.brep_optimized,0,'Optimized Geometry');
                fem.printElascticityResults();
            end
        end

        function fem = createSolver(obj, params)
            % createSolver
            % Given a parameter set, create and return the FEA solver object
            brep = obj.m_brepHandle(params);
            fem  = obj.m_solverHandle(brep);
        end

        %% OBJECTIVE & CONSTRAINTS
        function fx = evaluateObjective(obj, x, reset)
            % evaluateObjective
            % A wrapper for the objective function (e.g. compliance).
            % This function also keeps track of the parameter and objective
            % history in persistent variables, which are then exported to
            % the workspace for further analysis/plotting.

            persistent paramsHistory objectiveHistory feasibleHistory

            % Reset persistent variables if requested
            if nargin > 2 && reset
                paramsHistory = [];
                objectiveHistory = [];
                feasibleHistory = [];
                fx = [];
                return;
            end

            % Scale the parameters using the nominal values
            params = x .* obj.m_param0.value;

            % Evaluate the compliance objective relative to the initial baseline
            if strcmp(obj.m_objective,'compliance')
                fem = obj.solve(params);
                cx  = fem.computeCompliance();
                fx  = cx / obj.m_cx0;
            else
                disp(['Objective ' obj.m_objective ' is not implemented!']);
                fx = 0; % fallback
                cx = 0;
            end

            % Store history for analysis
            paramsHistory    = [paramsHistory;    params];
            objectiveHistory = [objectiveHistory; cx];
            % Mark feasibility if within bounds
            inBounds = all(params <= obj.m_param0.ub & params >= obj.m_param0.lb);
            feasibleHistory = [feasibleHistory; double(inBounds)];

            % Export to workspace
            assignin('base', 'paramsHistory', paramsHistory);
            assignin('base', 'objectiveHistory', objectiveHistory);
            assignin('base', 'feasibleHistory', feasibleHistory);

            % Optionally save .gif
            if (obj.m_exportGIF == 1)
                export_gifs();
            end
        end

        function [cineq, ceq] = evaluateConstraint(obj, x, reset)
            % evaluateConstraint
            % Dispatch to evaluate either area or perimeter constraint.

            if nargin < 3
                reset = 0;
            end

            if isfield(obj.m_constraints, 'area')
                [cineq, ceq] = obj.evaluateAreaConstraint(x, reset);
            elseif isfield(obj.m_constraints, 'perimeter')
                [cineq, ceq] = obj.evaluatePerimeterConstraint(x, reset);
            else
                disp('No constraint imposed!')
                cineq = [];
                ceq   = [];
            end
        end

        function [cineq, ceq] = evaluateAreaConstraint(obj, x, reset)
            % evaluateAreaConstraint
            % Evaluate the (in)equality constraint on the area.
            persistent areaHistory

            % Reset persistent variables if requested
            if nargin > 2 && reset
                areaHistory = [];
                cineq = [];
                ceq   = [];
                return;
            end

            % Scale the parameters
            params = x .* obj.m_param0.value;
            fem    = obj.createSolver(params);
            area   = fem.getArea();

            % If the user specified an inequality constraint
            if strcmp(obj.m_constraints.type, 'ineq')
                % cineq <= 0  =>  area <= target
                cineq = area / obj.m_targetArea - 1.0;
                ceq   = [];
            else
                % ceq = 0  =>  area == target
                cineq = [];
                ceq   = area / obj.m_targetArea - 1.0;
            end

            % Keep track of areas
            areaHistory = [areaHistory; area];
            assignin('base', 'areaHistory', areaHistory);
        end

        function [cineq, ceq] = evaluatePerimeterConstraint(obj, x, reset)
            % evaluatePerimeterConstraint
            % Evaluate the (in)equality constraint on the perimeter.
            persistent perimHistory

            % Reset persistent variables if requested
            if nargin > 2 && reset
                perimHistory = [];
                cineq = [];
                ceq   = [];
                return;
            end

            % Scale the parameters
            params = x .* obj.m_param0.value;
            fem    = obj.createSolver(params);
            perim  = fem.getPerimeter();

            % If the user specified an inequality constraint
            if strcmp(obj.m_constraints.type, 'ineq')
                % cineq <= 0  =>  perimeter <= target
                cineq = perim / obj.m_targetPerimeter - 1.0;
                ceq   = [];
            else
                % ceq = 0  =>  perimeter == target
                cineq = [];
                ceq   = perim / obj.m_targetPerimeter - 1.0;
            end

            % Keep track of perimeters
            perimHistory = [perimHistory; perim];
            assignin('base', 'perimHistory', perimHistory);
        end

        %% OPTIMIZATION ROUTINES
        function obj = optimize(obj)
            % optimize
            % Runs the optimization (fmincon, GlobalSearch, etc.) based on the chosen method.
            % Also updates the unified m_results structure with final info.

            scenarioId = 1;
            % Initialize structure to store solutions from certain solvers
            obj.m_feasibleExploredSolutions = struct('X',[],'Fval',[]);

            % Lower and upper bounds are relative factors (x = param/param0.value)
            LB = (obj.m_param0.lb) ./ obj.m_param0.value;
            UB = (obj.m_param0.ub) ./ obj.m_param0.value;

            % Starting guess in scaled space
            x0 = ones(size(LB));

            if obj.m_method == "FD"
                %===========================
                % Finite Differences (fmincon)
                %===========================
                opt = optimoptions('fmincon','Display','iter',...
                    'TolX',obj.m_terminationTolerance,...
                    'TolFun',obj.m_terminationTolerance,...
                    'ConstraintTolerance',obj.m_terminationTolerance,...
                    'FiniteDifferenceStepSize',obj.m_finiteDifferenceStepSize);

                [xMin, ~, flag, output] = fmincon(@obj.evaluateObjective, x0, [],[],...
                    [],[], LB,UB, @obj.evaluateConstraint, opt);

            elseif obj.m_method == "GS"
                %===========================
                % Global Search
                %===========================
                opt = optimoptions('fmincon',...
                    'TolX',obj.m_terminationTolerance,...
                    'TolFun',obj.m_terminationTolerance,...
                    'ConstraintTolerance',obj.m_terminationTolerance,...
                    'FiniteDifferenceStepSize',obj.m_finiteDifferenceStepSize);

                problem = createOptimProblem('fmincon','objective', @obj.evaluateObjective, ...
                    'x0', x0, 'lb', LB,'ub', UB, 'nonlcon', @obj.evaluateConstraint, 'options', opt);

                gs = GlobalSearch('MaxTime', 300, 'NumTrialPoints', 300, 'NumStageOnePoints', 100);
                [xMin, ~, flag, output, solutions] = run(gs, problem);
                obj.m_feasibleExploredSolutions = solutions;

            elseif obj.m_method == "MS"
                if isempty(obj.m_numLocalMultiStartProblems) || isnan(obj.m_numLocalMultiStartProblems)
                    disp('WARNING: number of local multi-start problems not set, setting default to 5.');
                    obj.m_numLocalMultiStartProblems = 5;
                end
                %===========================
                % Multi-Start
                %===========================
                opt = optimoptions('fmincon',...
                    'TolX',obj.m_terminationTolerance,...
                    'TolFun',obj.m_terminationTolerance,...
                    'ConstraintTolerance',obj.m_terminationTolerance,...
                    'FiniteDifferenceStepSize',obj.m_finiteDifferenceStepSize);

                problem = createOptimProblem('fmincon','objective', @obj.evaluateObjective, ...
                    'x0', x0, 'lb', LB,'ub', UB, 'nonlcon', @obj.evaluateConstraint, 'options', opt);

                ms = MultiStart('UseParallel',0);
                [xMin, ~, flag, output, solutions] = run(ms, problem, obj.m_numLocalMultiStartProblems);
                obj.m_feasibleExploredSolutions = solutions;

            elseif obj.m_method == "RS"
                %===========================
                % Random Search
                %===========================
                if isempty(obj.m_numRandomSamples) || isnan(obj.m_numRandomSamples)
                    disp('WARNING: number of samples not set, setting default to 10.');
                    obj.m_numRandomSamples = 10;
                end
                if isempty(obj.m_numRandomFuncCounts) || isnan(obj.m_numRandomFuncCounts)
                    disp('WARNING: number of max function evaluations is not set, setting default to 100.');
                    obj.m_numRandomFuncCounts = 100;
                end

                nSamples = obj.m_numRandomSamples;
                nFuncCount = obj.m_numRandomFuncCounts;
                bestObj  = Inf;
                bestX    = [];

                solutionsRS = [];  % or pre-allocate a struct array
                % Track number of objective evaluations
                funcCount = 0;
                feasSolCount = 0;
                while feasSolCount < nSamples && funcCount < nFuncCount
                    % Random point in [LB, UB]
                    xCandidate = LB + rand(size(LB)).*(UB - LB);

                    % Evaluate constraints
                    [c, ceq] = obj.evaluateConstraint(xCandidate);

                    % Evaluate objective
                    fval = obj.evaluateObjective(xCandidate);
                    funcCount = funcCount + 1;

                    % Check feasibility
                    isFeasible = all(c <= 0) && all(abs(ceq) <= 1e-6);
                    if isFeasible
                        feasSolCount = feasSolCount + 1;
                        solutionsRS(feasSolCount).X    = xCandidate; %#ok        % scaled param
                        solutionsRS(feasSolCount).Fval = fval;   %#ok            % objective
                        if (fval < bestObj)
                            bestObj = fval;
                            bestX   = xCandidate;
                        end
                    end
                end

                % Check if we found a feasible solution
                if ~isempty(bestX)
                    xMin = bestX;
                    flag = 1;  % success
                else
                    xMin = x0;
                    flag = -1; % no feasible solution found
                end

                % Create an output struct for consistency
                output.funcCount = funcCount;
                output.message   = 'Random Search completed';
                obj.m_feasibleExploredSolutions = solutionsRS;
            else
                %===========================
                % Fallback to Finite Differences
                %===========================
                disp("Unknown method, setting to FD as default.")
                opt = optimoptions('fmincon','Display','iter',...
                    'TolX',obj.m_terminationTolerance,...
                    'TolFun',obj.m_terminationTolerance,...
                    'ConstraintTolerance',obj.m_terminationTolerance,...
                    'FiniteDifferenceStepSize',obj.m_finiteDifferenceStepSize);

                [xMin, ~, flag, output] = fmincon(@obj.evaluateObjective, x0, [],[],...
                    [],[], LB,UB, @obj.evaluateConstraint, opt);
            end

            %===========================
            % Store results
            %===========================
            obj.m_results.flag = flag;

            % Record initial baseline
            obj.m_results.initialParams     = obj.m_param0.value;
            obj.m_results.initialArea       = obj.m_area0;
            obj.m_results.initialPerimeter  = obj.m_perim0;
            obj.m_results.initialCompliance = obj.m_cx0;
            obj.m_results.initialDef        = obj.m_maxDef0;
            obj.m_results.initialStress     = obj.m_maxStress0;

            % Record final solution
            obj.m_optimalParams = xMin .* obj.m_param0.value;
            fem = obj.solve(obj.m_optimalParams);

            obj.m_results.finalParams     = obj.m_optimalParams;
            obj.m_results.finalArea       = fem.getArea();
            obj.m_results.finalPerimeter  = fem.getPerimeter();
            obj.m_results.finalCompliance = fem.computeCompliance();
            obj.m_results.finalDef        = fem.m_maxDef(scenarioId);
            obj.m_results.finalStress     = fem.m_maxStress(scenarioId);

            % Approx # of FEAs (fmincon usually calls objective + grad)
            if isfield(output,'funcCount')
                % For random search, we used our own counting
                feCount = output.funcCount;
            else
                % For fmincon-based methods
                feCount = output.funcCount * (1 + length(obj.m_param0.value));
            end
            obj.m_results.nFEAs = feCount;

            % Store final solver
            obj.m_solverFinal = fem;

            % Harvest iteration history from the workspace
            obj.m_results.history.x        = evalin('base', 'paramsHistory');
            obj.m_results.history.f        = evalin('base', 'objectiveHistory');
            obj.m_results.history.feasible = evalin('base', 'feasibleHistory');

            if evalin('base', 'exist(''areaHistory'', ''var'')')
                obj.m_results.history.area = evalin('base', 'areaHistory');
            end
            if evalin('base', 'exist(''perimHistory'', ''var'')')
                obj.m_results.history.perim = evalin('base', 'perimHistory');
            end

            if ~obj.m_testMode
                % Display basic summary
                disp('Exit Flag: '); disp(obj.m_results.flag);
                disp('Initial Params: '); disp(obj.m_results.initialParams);
                disp('Final Params: ');   disp(obj.m_results.finalParams);

                disp('Initial/Final Area: ');
                disp([obj.m_results.initialArea obj.m_results.finalArea]);

                disp('Initial/Final Perimeter: ');
                disp([obj.m_results.initialPerimeter obj.m_results.finalPerimeter]);

                disp('Initial/Final Compliance: ');
                disp([obj.m_results.initialCompliance obj.m_results.finalCompliance]);

                disp('Initial/Final Def: ');
                disp([obj.m_results.initialDef obj.m_results.finalDef]);

                disp('Initial/Final Stress: ');
                disp([obj.m_results.initialStress obj.m_results.finalStress]);

                disp('#FEAs:'); disp(obj.m_results.nFEAs);

                % Simple post-processing plot(s)
                if obj.m_method == "FD"
                    % For FD, let's assume we want to plot typical "Convergence"
                    obj = obj.plotConvergence();
                else
                    obj = obj.plotParetoSpace();
                end
            end
        end

        function obj = plotConvergence(obj)
            % plotConvergence
            % Plots objective and constraint history over iterations.
            plt = PlotId;
            figure(plt.convergence);
            set(gcf, 'Name', 'Convergence');

            nIter = numel(obj.m_results.history.f);

            % Identify whether we are using area or perimeter
            if isfield(obj.m_constraints, 'area')
                consVals = obj.m_results.history.area(1:nIter);
                consLabel = 'Area ($A$)';
            elseif isfield(obj.m_constraints, 'perimeter')
                consVals = obj.m_results.history.perim(1:nIter);
                consLabel = 'Perimeter ($P$)';
            else
                disp('No constraint imposed!');
                consVals = zeros(1,nIter);
                consLabel = 'No Constraint';
            end

            % Plot objective
            f = obj.m_results.history.f(1:nIter);
            plot(1:nIter, f, '-b','LineWidth',2);
            hold on; grid on;
            xlabel('Iteration');

            if strcmp(obj.m_objective,'compliance')
                ylabel('Compliance $C$');
            else
                ylabel('Objective $\varphi$');
            end
            axis tight;

            % Plot constraint with a secondary y-axis
            yyaxis right
            plot(1:nIter, consVals, '--r','LineWidth',2);
            ylabel(consLabel);

            hold off;
        end

        function obj = plotParetoSpace(obj)
            % plotParetoSpace
            % Scatter plots the explored solutions in terms of area (or perimeter)
            % vs. objective (e.g. compliance ratio). Solutions come from the
            % GlobalSearch (or other multi-start) results. The optimum point
            % is highlighted in red.

            plt = PlotId;

            % Check whether 'area' or 'perimeter' is constrained
            useArea = false;
            usePerimeter = false;
            if isfield(obj.m_constraints, 'area')
                useArea = true;
            elseif isfield(obj.m_constraints, 'perimeter')
                usePerimeter = true;
            else
                warning('No constraint on area or perimeter found!');
            end

            % Number of solutions found (local minima) by GlobalSearch
            solutions = obj.m_feasibleExploredSolutions;
            nSolutions = numel(solutions);

            if nSolutions == 0
                warning('No solutions recorded in m_feasibleExploredSolutions!');
                return;
            end

            % Prepare arrays for x (constraint) and y (objective)
            xVals = zeros(nSolutions,1);
            yVals = zeros(nSolutions,1);

            % Loop through each solution, re-solve to get area/perimeter
            for i = 1:nSolutions
                % Each solution has X (scaled param set) and Fval (objective)
                xCandidate = solutions(i).X;
                fCandidate = solutions(i).Fval;  % This should be compliance ratio if that’s your objective

                % Convert scaled param to physical param
                paramCandidate = xCandidate .* obj.m_param0.value;

                % Create and solve the FEA for that param set
                femCandidate = obj.createSolver(paramCandidate);
                femCandidate = femCandidate.preProcess();
                femCandidate = femCandidate.solve();
                femCandidate = femCandidate.postProcess();

                % Extract the constraint measure (area or perimeter)
                if useArea
                    xVals(i) = femCandidate.getArea();
                elseif usePerimeter
                    xVals(i) = femCandidate.getPerimeter();
                else
                    xVals(i) = NaN;  % fallback
                end

                % Objective is stored in .Fval
                yVals(i) = fCandidate*obj.m_results.initialCompliance;
            end

            % Now, get the final (optimal) point’s area/perimeter and objective ratio
            if useArea
                xOpt = obj.m_results.finalArea;
            elseif usePerimeter
                xOpt = obj.m_results.finalPerimeter;
            else
                xOpt = NaN;
            end

            % If objective is compliance ratio, we can compute it from final results:
            % ratio = finalCompliance
            yOpt = obj.m_results.finalCompliance;

            % Create figure and scatter plot
            figure(plt.pareto_front);
            set(gcf, 'Name', 'Pareto Space');
            hold on; grid on;
            scatter(xVals, yVals, 40, 'b', 'filled', ...
                'DisplayName','Explored Solutions');
            scatter(xOpt, yOpt, 80, 'r', 'filled', ...
                'DisplayName','Optimum');

            % Labels
            if useArea
                xlabel('Area');
            elseif usePerimeter
                xlabel('Perimeter');
            else
                xlabel('Constraint Value');
            end

            ylabel('Objective (Compliance)');
            legend('Location','best');

            hold off;
        end


    end
end
