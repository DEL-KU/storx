%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discription:                                                              %
% This is an FEA code for 2D thermoelasticity problems. It inherits from the%
% simulation2d class and provides methods for pre-processing, solving, and  %
% post-processing the simulation. Thermal and elasticity solvers are        %
% utilized to solve weakly coupled thermoelasticity problems where we first %
% solve the thermal problem to obtain the temperature field,                %
% then use the temperature field to compute thermal strain and stress, and  %
% finally solve the elasticity problem to obtain the displacement field.    %
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

classdef fea2d_thermoelasticity < simulation2d
    properties(GetAccess = 'public', SetAccess = 'public')
        m_thermalStrain; % thermal strain (nely x nelx x 8)
        m_thermalStress; % thermal stress (nely x nelx)
        m_TReference = 0; % Reference Temperature Field
        m_fThermal; % thermal load based on thermal stress
        m_dThermalStressdT; % gradient of thermal stress with respect to temperature
        m_thermalSolver; % thermal conduction solver
        m_elasticitySolver; % elasticity solver
        m_constantDeltaT = 1;  % constant delta T in the entire domain
        m_isConstantDeltaT = false; % flag for constant delta T
    end
    methods
        function obj = fea2d_thermoelasticity(thermalSolver, elasticitySolver,TReference,...
                brep,numElements,materials,vectorize,...
                numScenarios,interpolation,uniformGrid)
            % set default values
            if nargin < 7, vectorize = false;end
            if nargin < 8, numScenarios = 1;end
            if nargin < 9, interpolation = 'none';end
            if nargin < 10, uniformGrid = 0;end

            % check inputs
            if ~isa(thermalSolver, 'fea2d_thermal')
                error('thermalSolver must be an instance of fea2d_thermal class!');
            end
            if ~isa(elasticitySolver, 'fea2d_elasticity')
                error('elasticitySolver must be an instance of fea2d_elasticity class!');
            end

            % construct the object
            obj = obj@simulation2d(brep,numElements,materials,...
                interpolation,numScenarios,uniformGrid); % call superclass

            obj.m_vectorize = vectorize;
            obj.m_thermalSolver = thermalSolver;
            obj.m_elasticitySolver = elasticitySolver;

            obj.m_TReference = TReference; % reference temperature

            obj.m_fThermal = zeros(size(obj.m_elasticitySolver.getForce()));
            obj.m_thermalStrain = zeros(obj.m_elasticitySolver.m_ny,...
                obj.m_elasticitySolver.m_nx,8);
            obj.m_thermalStress = zeros(obj.m_elasticitySolver.m_ny,...
                obj.m_elasticitySolver.m_nx,8);
        end

        function obj = preProcess(obj)
            % pre-process the thermal and elasticity solvers
            obj.m_thermalSolver = obj.m_thermalSolver.preProcess();
            obj.m_elasticitySolver = obj.m_elasticitySolver.preProcess();

        end

        function obj = postProcess(obj)
            % post-process the thermal and elasticity solvers
            obj.m_thermalSolver = obj.m_thermalSolver.postProcess();
            obj.m_elasticitySolver = obj.m_elasticitySolver.postProcess();
        end


        function obj = setConstantDeltaT(obj,deltaT)
            % set a constant delta T for the entire domain
            obj.m_constantDeltaT = deltaT;
            obj.m_isConstantDeltaT = true;
        end
        %% PENALTY CONTINUATION
        function obj = performPenaltyContinuation(obj)
            % perform penalty continuation for both thermal and elasticity solvers
            obj.m_thermalSolver = obj.m_thermalSolver.performPenaltyContinuation();
            obj.m_elasticitySolver = obj.m_elasticitySolver.performPenaltyContinuation();
        end

        function obj = setPenaltyFactor(obj,penalty)
            % set penalty factor for both thermal and elasticity solvers
            obj.m_thermalSolver = obj.m_thermalSolver.setPenaltyFactor(penalty);
            obj.m_elasticitySolver = obj.m_elasticitySolver.setPenaltyFactor(penalty);
        end

        function obj = solve(obj)
            % update designs
            obj = obj.updateDesign();
            % solve thermal
            obj.m_thermalSolver = obj.m_thermalSolver.solve();
            obj.m_thermalSolver = obj.m_thermalSolver.postProcess();
            % % compute structural loads from thermal expansion
            obj = obj.computeThermalForce();
            % % update structural loads
            fStructural = obj.m_elasticitySolver.getForce();
            fThermoelastic =  fStructural + obj.m_fThermal;
            obj.m_elasticitySolver = obj.m_elasticitySolver.setForce(fThermoelastic);
            % solve elasticity
            obj.m_elasticitySolver = obj.m_elasticitySolver.solve();
            % reset the structural load
            obj.m_elasticitySolver = obj.m_elasticitySolver.setForce(fStructural);
        end

        function obj = computeThermalForce(obj)
            if obj.m_vectorize
                obj = obj.computeThermalForce_vectorized();
            else
                obj = obj.computeThermalForce_loop();
            end
        end

        function obj = computeThermalForce_vectorized(obj)
            % compute thermal force using vectorized operations
            % this is the vectorized version of the thermal force computation
            % it is faster than the loop version and is suitable for large problems
            % initialize thermal strain and stress

            nelx = obj.m_elasticitySolver.m_nx;
            nely = obj.m_elasticitySolver.m_ny;
            nDOF = obj.m_elasticitySolver.m_numDOFs;
            numScenarios = obj.m_elasticitySolver.m_numScenarios;
            obj.m_thermalStrain = zeros(nely, nelx, 8);
            obj.m_fThermal = sparse(nDOF, numScenarios); % initialize thermal stress load
            obj.m_dThermalStressdT = zeros(nely, nelx, 8); % initialize gradient of thermal stress with respect to temperature
            % Constants for nodal positions
            xNodes = [0, 1, 1, 0] * obj.m_elasticitySolver.m_hx;
            yNodes = [0, 0, 1, 1] * obj.m_elasticitySolver.m_hy;
            xi = 0; eta = 0;
            J = obj.m_elasticitySolver.Jacobian(xNodes, yNodes, xi, eta);
            [~, gradN] = obj.m_elasticitySolver.QuadShapeFunction(xi, eta);
            B = J' \ gradN;
            Z = zeros(1, 4);
            Bmat = [B(1, :) Z; Z B(2, :); B(2, :) B(1, :)];
            Phi = [1 1 0];


            ve = obj.m_elasticitySolver.m_ve; % element volume
            TReference = obj.m_TReference; % reference temperature

            % Vectorize over scenarios
            for scenarioId = 1:numScenarios
                TemperatureField = obj.m_thermalSolver.getTemperature(scenarioId);

                % Logical indexing to select valid elements
                validElems = obj.m_elasticitySolver.m_existingElems;
                [elyGrid, elxGrid] = find(validElems);

                % Get material indices for all elements
                matIdGrid = obj.m_elasticitySolver.m_materialIndices(sub2ind(size(validElems), elyGrid, elxGrid));

                % Node indices for all elements
                nodes = [(elxGrid - 1) * (nely + 1) + elyGrid, elxGrid * (nely + 1) + elyGrid, ...
                    elxGrid * (nely + 1) + elyGrid + 1, (elxGrid - 1) * (nely + 1) + elyGrid + 1];

                % Compute deltaT for each element
                if obj.m_isConstantDeltaT
                    deltaTGrid = obj.m_constantDeltaT * ones(size(elxGrid));
                else
                    TelemGrid = 0.25 * (TemperatureField(sub2ind(size(TemperatureField), elyGrid, elxGrid)) + ...
                        TemperatureField(sub2ind(size(TemperatureField), elyGrid, elxGrid + 1)) + ...
                        TemperatureField(sub2ind(size(TemperatureField), elyGrid + 1, elxGrid + 1)) + ...
                        TemperatureField(sub2ind(size(TemperatureField), elyGrid + 1, elxGrid)));
                    deltaTGrid = TelemGrid - TReference;
                end

                % Calculate thermal strain and stress for all elements
                for idx = 1:length(elxGrid)
                    matId = matIdGrid(idx);
                    alpha = obj.m_elasticitySolver.m_materials(matId).alpha;
                    deltaTelem = deltaTGrid(idx);
                    thermalStrain = alpha * deltaTelem * Phi';  % thermal strain

                    D = obj.m_elasticitySolver.getElasticityMatrix(matId); % elasticity matrix

                    thermalStress = Bmat' * D * thermalStrain * ve; % thermal stress

                    x = obj.m_design(elyGrid(idx), elxGrid(idx));   % design
                    interpCoeff = obj.m_thermalSolver.getInterpolationCoefficient(x); % interpolation coefficient

                    order = reshape([1:4; 5:8], 1, 8); % alternate u and v
                    thermalStress = thermalStress(order);

                    % Store thermal stress in element
                    obj.m_thermalStress(elyGrid(idx), elxGrid(idx), :) = interpCoeff * thermalStress;

                    % gradient of thermal stress with respect to temperature
                    obj.m_dThermalStressdT(elyGrid(idx), elxGrid(idx), :)  = 0.25*obj.m_thermalStress(elyGrid(idx), elxGrid(idx), :) /max(deltaTelem,1e-9);

                    % Compute DOF for each element
                    dof = [2 * nodes(idx, :) - 1; 2 * nodes(idx, :)];
                    dof = reshape(dof, 1, 8);

                    % Accumulate thermal force
                    obj.m_fThermal(dof, scenarioId) = obj.m_fThermal(dof, scenarioId) + thermalStress;
                end
            end
        end


        function obj = computeThermalForce_loop(obj)
            % compute thermal force using loops
            % this is the non-vectorized version of the thermal force computation
            % it is slower but can be useful for debugging or small problems
            % where vectorization is not necessary
            % initialize thermal strain and stress

            nelx = obj.m_elasticitySolver.m_nx;
            nely = obj.m_elasticitySolver.m_ny;
            nDOF = obj.m_elasticitySolver.m_numDOFs;
            obj.m_thermalStrain = zeros(obj.m_elasticitySolver.m_ny,...
                obj.m_elasticitySolver.m_ny,8);
            obj.m_fThermal = sparse(nDOF ,obj.m_elasticitySolver.m_numScenarios); % initialize thermal stress load
            obj.m_dThermalStressdT = zeros(obj.m_elasticitySolver.m_ny, obj.m_elasticitySolver.m_ny, 8); % initialize gradient of thermal stress with respect to temperature
            xNodes = [0,1,1,0]*obj.m_elasticitySolver.m_hx;
            yNodes = [0,0,1,1]*obj.m_elasticitySolver.m_hy;
            xi = 0;
            eta = 0;
            J = obj.m_elasticitySolver.Jacobian(xNodes,yNodes,xi,eta);
            [~,gradN] = obj.m_elasticitySolver.QuadShapeFunction(xi,eta);
            B = J'\gradN;
            Z = zeros(1,4);
            Bmat = [B(1,:) Z; Z B(2,:); B(2,:) B(1,:)];
            Phi = [1 1 0];
            for scenarioId = 1:obj.m_elasticitySolver.m_numScenarios
                TemperatureField = obj.m_thermalSolver.getTemperature(scenarioId);
                for elx = 1:nelx
                    for ely = 1:nely
                        if (~obj.m_elasticitySolver.m_existingElems(ely,elx)),continue;end

                        matId = obj.m_elasticitySolver.m_materialIndices(ely,elx);
                        nodes = [((elx-1)*(nely+1) + ely),...
                            ((elx)*(nely+1) + ely),...
                            ((elx)*(nely+1) + ely+1),...
                            ((elx-1)*(nely+1) + ely+1)];

                        if (obj.m_isConstantDeltaT)
                            deltaTelem = obj.m_constantDeltaT;
                        else
                            % average nodal temperatures at each element
                            Telem = 0.25*(TemperatureField(ely, elx)+...
                                TemperatureField(ely, elx+1) + ...
                                TemperatureField(ely+1, elx+1) + ...
                                TemperatureField(ely+1, elx));

                            deltaTelem = Telem - obj.m_TReference;
                        end

                        alpha = obj.m_elasticitySolver.m_materials(matId).alpha;
                        thermalStrain = alpha*deltaTelem*Phi';

                        D = obj.m_elasticitySolver.getElasticityMatrix(matId);

                        thermalStress= Bmat'*D*thermalStrain ...
                            *obj.m_elasticitySolver.m_ve;

                        x = obj.m_design(ely,elx);
                        interpCoeff = obj.m_thermalSolver.getInterpolationCoefficient(x);

                        order = reshape([1:4;5:8],1,8); % alternate u and v
                        thermalStress = thermalStress(order);
                        obj.m_thermalStress(ely,elx,:) = interpCoeff*thermalStress;

                        % gradient of thermal stress with respect to temperature
                        obj.m_dThermalStressdT(ely,elx,:) = 0.25*obj.m_thermalStress(ely,elx,:)/deltaTelem;

                        dof = [2*nodes-1; 2*nodes];
                        dof = reshape(dof,1,8);
                        obj.m_fThermal(dof,scenarioId) = obj.m_fThermal(dof,scenarioId) + thermalStress;
                    end
                end
            end
        end

        %% DESIGN UPDATE
        function obj = updateDesign(obj)
            % update the design in both thermal and elasticity solvers
            % this is called before solving the problem
            obj.m_thermalSolver = obj.m_thermalSolver.setDesign(obj.m_design);
            obj.m_elasticitySolver = obj.m_elasticitySolver.setDesign(obj.m_design);
        end

        %% PLOTTING
        function obj = plotBoundaryCondition(obj)
            plt = PlotId;
            plbc = PlotBC;
            for scenarioId = 1:obj.m_numScenarios
                legend_fields = {};
                legend_labels = {};
                % plot grid mesh, sets outside value to NaN, so they are ignored in the plot
                obj.plotGeometry(plt.loading+scenarioId,0);
                set(gcf, 'Name', strjoin({'Boundary Condition',num2str(scenarioId)},' '));hold on;

                %% Thermal
                % mark all fixed nodes
                numNodes = (obj.m_nx+1)*(obj.m_ny+1);
                X = zeros(1,numNodes); Y = zeros(1,numNodes);
                hx = obj.m_boxSizes(1)/obj.m_nx;
                hy = obj.m_boxSizes(2)/obj.m_ny;
                for nodex = 1:obj.m_nx+1
                    for nodey = 1:obj.m_ny+1
                        nodeId = ((nodex-1)*(obj.m_ny+1) + nodey);
                        X(nodeId) = obj.m_boundingBox(1,1) + (nodex-1)*hx;
                        Y(nodeId) = obj.m_boundingBox(2,1) + (nodey-1)*hy;
                    end
                end
                X_TFixed = X(obj.m_thermalSolver.m_fixedDOFs);
                Y_TFixed = Y(obj.m_thermalSolver.m_fixedDOFs);
                if ~isempty(obj.m_thermalSolver.m_fixedDOFs)
                    fixed_T = plot(X_TFixed',Y_TFixed',plbc.fixed_T.marker,'MarkerEdgeColor', plbc.fixed_T.color,'MarkerSize',10);
                    hold on;
                    % legend
                    legend_fields = [legend_fields;fixed_T]; %#ok
                    legend_labels = [legend_labels,'fixed $T$']; %#ok
                end

                % mark all flux nodes
                nodes = obj.m_thermalSolver.m_forcedNodes{scenarioId};
                if ~isempty(nodes)
                    nodes = nodes(1:2:end);
                    X_forced = X(nodes);
                    Y_forced = Y(nodes);
                    Fx = obj.m_thermalSolver.m_f(nodes,scenarioId);
                    Fy = obj.m_thermalSolver.m_f(nodes,scenarioId);
                    scale = 0.1*obj.m_modelScale;
                    normF = sqrt(Fx.^2 + Fy.^2);
                    for i = 1:numel(nodes)
                        seg = obj.m_thermalSolver.m_fluxSegs{scenarioId}(i);
                        normal = obj.m_thermalSolver.normalOfSegment(seg);
                        Fx(i) = normal(1)*Fx(i);
                        Fy(i) = normal(2)*Fy(i);
                    end
                    if (normF > 0)
                        Fx = Fx./normF;
                        Fy = Fy./normF;
                        Fx = Fx(:)';
                        Fy = Fy(:)';
                        plot(X_forced,Y_forced,'Marker','diamond','MarkerFaceColor',plbc.flux.color);hold on;
                        start = [X_forced; Y_forced];
                        stop = start + scale*[Fx; Fy];
                        obj.drawArrow(start',stop',plbc.flux.color,2);hold on;
                        % legend
                        flux = plot(NaN, NaN,plbc.flux.marker,'MarkerEdgeColor',plbc.flux.color,'MarkerFaceColor', plbc.flux.color);
                        legend_fields = [legend_fields;flux ]; %#ok
                        legend_labels = [legend_labels,'heat flux']; %#ok
                    end
                end
                %%
                % internal heat
                ih = obj.m_thermalSolver.m_internalHeat(scenarioId);
                if abs(ih) > 0
                    % Get axis limits
                    x_limits = xlim;
                    y_limits = ylim;
                    start = obj.findPointinEmptyRegion(x_limits,y_limits,[X_TFixed;Y_TFixed]);
                    scale = 0.1*obj.m_modelScale;
                    nPoints = 8;
                    % Generate angles from 0 to 2*pi
                    theta = linspace(0, 2*pi, nPoints+1);
                    theta(end) = [];
                    % Parametric equations for the circle
                    stop = start + scale*[cos(theta);sin(theta)];
                    start = start + 0.2*scale*[cos(theta);sin(theta)];
                    if ih < 0
                        tmp = stop;
                        stop = start;
                        start = tmp;
                    end
                    obj.drawArrow(start',stop',plbc.internal_heat.color,4);hold on;
                    % legend
                    internal_heat = plot(NaN, NaN,plbc.internal_heat.marker,'MarkerEdgeColor',plbc.internal_heat.color,'MarkerFaceColor', plbc.internal_heat.color);
                    legend_fields = [legend_fields;internal_heat ]; %#ok
                    legend_labels = [legend_labels,'internal heat']; %#ok
                end

                %% Elasticity
                % mark all fixed nodes
                numNodes = (obj.m_nx+1)*(obj.m_ny+1);
                X = zeros(1,numNodes); Y = zeros(1,numNodes);
                hx = obj.m_boxSizes(1)/obj.m_nx;
                hy = obj.m_boxSizes(2)/obj.m_ny;
                for nodex = 1:obj.m_nx+1
                    for nodey = 1:obj.m_ny+1
                        nodeId = ((nodex-1)*(obj.m_ny+1) + nodey);
                        X(nodeId) = obj.m_boundingBox(1,1) + (nodex-1)*hx;
                        Y(nodeId) = obj.m_boundingBox(2,1) + (nodey-1)*hy;
                    end
                end
                index = (rem(obj.m_elasticitySolver.m_fixedDOFs,2)==1);
                fixedXNodes = (obj.m_elasticitySolver.m_fixedDOFs(index)-1)/2+1;
                X_xFixed = X(fixedXNodes);
                Y_xFixed = Y(fixedXNodes);
                if ~isempty(fixedXNodes)
                    fixed_U = plot(X_xFixed',Y_xFixed',plbc.fixed_U.marker,'MarkerEdgeColor', plbc.fixed_U.color,'LineWidth',1);
                    hold on;
                    % legend
                    legend_fields = [legend_fields;fixed_U ]; %#ok
                    legend_labels = [legend_labels,'fixed $u$']; %#ok
                end
                index = (rem(obj.m_elasticitySolver.m_fixedDOFs,2)==0);
                fixedYNodes = (obj.m_elasticitySolver.m_fixedDOFs(index))/2;
                X_yFixed = X(fixedYNodes);
                Y_yFixed = Y(fixedYNodes);
                if ~isempty(fixedYNodes)
                    fixed_V = plot(X_yFixed',Y_yFixed',plbc.fixed_V.marker,'MarkerEdgeColor', plbc.fixed_V.color,'LineWidth',1);
                    hold on;
                    % legend
                    legend_fields = [legend_fields;fixed_V ]; %#ok
                    legend_labels = [legend_labels,'fixed $v$']; %#ok
                end

                % mark all forced nodes
                nodes = obj.m_elasticitySolver.m_forcedNodes{scenarioId};
                if ~isempty(nodes)
                    nodes = nodes(1:2:end);
                    X_forced = X(nodes);
                    Y_forced = Y(nodes);
                    Fx = obj.m_elasticitySolver.m_f(2*nodes-1,scenarioId);
                    Fy = obj.m_elasticitySolver.m_f(2*nodes,scenarioId);
                    scale = 0.1*obj.m_modelScale;
                    normF = sqrt(Fx.^2 + Fy.^2);
                    if (normF > 0)
                        Fx = Fx./normF;
                        Fy = Fy./normF;
                        Fx = Fx(:)';
                        Fy = Fy(:)';
                        plot(X_forced,Y_forced,'Marker','diamond','MarkerFaceColor',plbc.force.color);hold on;
                        start = [X_forced; Y_forced];
                        stop = start + scale*[Fx; Fy];
                        obj.drawArrow(start',stop',plbc.force.color,2);hold on;
                        % legend
                        force = plot(NaN, NaN,plbc.force.marker,'MarkerEdgeColor',plbc.force.color,'MarkerFaceColor', plbc.force.color);
                        legend_fields = [legend_fields;force ]; %#ok
                        legend_labels = [legend_labels,'force']; %#ok
                    end
                end
                % acceleration
                bx = obj.m_elasticitySolver.m_acceleration(scenarioId,1);
                by = obj.m_elasticitySolver.m_acceleration(scenarioId,2);
                bNorm = sqrt(bx^2+by^2);
                if bNorm > 1e-6
                    % Get axis limits
                    x_limits = xlim;
                    y_limits = ylim;
                    scale = 0.15*obj.m_modelScale;
                    start = obj.findPointinEmptyRegion(x_limits,y_limits,[[X_xFixed,X_yFixed,X_TFixed];[Y_xFixed,Y_yFixed,Y_TFixed]]);
                    stop = start + scale*[bx/bNorm; by/bNorm];
                    obj.drawArrow(start',stop',plbc.acceleration.color,4);hold on;
                    % legend
                    acceleration = plot(NaN, NaN,plbc.acceleration.marker,'MarkerEdgeColor',plbc.acceleration.color,'MarkerFaceColor', plbc.acceleration.color);
                    legend_fields = [legend_fields;acceleration ]; %#ok
                    legend_labels = [legend_labels,'body force']; %#ok
                end


                legend(legend_fields,legend_labels, ...
                    'Location', 'northeastoutside');
                pbaspect(obj.m_boxSizes);axis on;axis tight;
            end
        end
    end

end