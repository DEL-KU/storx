%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                              %
% This is the abstract class for 2D density-based topology optimization.    %
% It is used to define the basic properties and methods for density-based   %
% topology optimization problems.                                           %
%                                                                           %
% The class inherits from topopt2d, and is designed to be inherited by      %
% specific density-based                                                    %
% topology optimization classes, which will implement the abstract methods. %
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

classdef (Abstract) density2d < topopt2d
    properties
        m_update; % update method (e.g., OC, MMA, GCMMA)
        m_ocOptimizer; % OC optimizer
        m_ocExponent = 2;   % exponent for OC
        m_mmaOptimizer; % MMA optimizer
        m_gcmmaOptimizer; % GCMMA optimizer
        m_targetVolFrac; % targt volume fraction
        m_stageDesign; % keep track of chain rule
        m_scale; % scaling factor for gradients
    end
    methods (Abstract)
        obj = saveHistory(obj)
    end
    methods
        function obj = density2d(solver,objective,constraints, mfgConstraints, ...
                update, maxNumIters,exportGIF,testMode)

            % construct
            obj = obj@topopt2d(solver,objective,constraints,mfgConstraints, ...
                maxNumIters,exportGIF,testMode); % call superclass
            for consId = 1 : obj.m_numConstraints
                if (isa(obj.m_constraints{consId},'volume'))
                    obj.m_targetVolFrac = obj.m_constraints{consId}.m_upperBound;
                    break;
                elseif (isa(obj.m_constraints{consId},'localVolume'))
                    obj.m_targetVolFrac = obj.m_constraints{consId}.m_upperBound;
                    break;
                elseif (isa(obj.m_constraints{consId},'activeVolume'))
                    obj.m_targetVolFrac = obj.m_constraints{consId}.m_upperBound;
                    break;
                end
            end

            obj.m_update = update;
            obj.m_x = obj.m_solver.m_existingElems*obj.m_targetVolFrac;
            obj.m_solver = obj.m_solver.setDesign(obj.m_x);
            

            for consId = 1 : obj.m_numMfgConstraints
                if (isa(obj.m_mfgConstraints{consId},'physicalDensity'))
                    if (strcmp(obj.m_update,'OC') == 1)
                        fprintf('\033[1m');
                        fprintf('=====================\n');
                        fprintf('WARNING: The physical density projection filter requires MMA or GCMMA optimizers!\n');
                        fprintf('WARNING: Changing the optimizer to MMA and maximum projection parameter set to 4!\n');
                        fprintf('=====================\n');
                        fprintf('\033[0m');
                        obj.m_update = 'MMA';
                        obj.m_mfgConstraints{consId} = obj.m_mfgConstraints{consId}.setParameters(1,0.5,4);
                    elseif (strcmp(obj.m_update,'MMA') == 1)
                        fprintf('\033[1m');
                        disp('=====================')
                        disp('WARNING: Setting the maximum projection parameter set to 4!');
                        disp('=====================')
                        fprintf('\033[0m');
                       obj.m_mfgConstraints{consId} = obj.m_mfgConstraints{consId}.setParameters(1,0.5,4);
                    end
                
                end
            end


            if (strcmp(obj.m_update,'OC') == 1)
                obj.m_ocOptimizer = ocOptimizer(0,1e9,0.2,0.001,1,obj.m_solver.m_numExistingElems*obj.m_targetVolFrac);
                obj.m_changeTarget = 0.02;
            elseif (strcmp(obj.m_update,'MMA') == 1)
                obj.m_mmaOptimizer = mmaOptimizer(obj.m_solver.m_nx*obj.m_solver.m_ny,obj.m_numConstraints,0,1000,1,obj.m_x(:));
                obj.m_changeTarget = 0.002;
            elseif (strcmp(obj.m_update,'GCMMA') == 1)
                obj.m_gcmmaOptimizer = gcmmaOptimizer(obj.m_solver.m_nx*obj.m_solver.m_ny,obj.m_numConstraints,0,1000,1,obj.m_x(:));
                obj.m_changeTarget = 0.002;
            end
        end

        function obj = setChangeTolerance(obj,targetTol)
            obj.m_changeTarget = targetTol;
        end

        function obj = optimize(obj)
            obj.m_flag = -1; % optimization flag, intialize to failed
            % Main optimization loop
            for iter = 1:obj.m_maxNumIters
                obj.m_iter = iter;
                %% filter density
                obj= obj.filterDensity();
                %% solve
                obj = obj.solve();
                %% evaluate and gradient
                obj = obj.evaluate();
                obj = obj.gradient();
                %% filter sensitivity
                obj = obj.filterSensitivity();
                %% history and display
                obj = obj.saveHistory();
                if (~obj.m_testMode)
                    obj = obj.printResults();  % print results if not in test mode
                    obj.plotPseudoDensity(); % plot design if not in test mode
                    if (obj.m_exportGIF), export_gifs();end
                end
                %% check convergence
                [obj,isConverged] = obj.checkConvergence();
                if (isConverged),obj.m_flag = 1; break; end
                %% update
                [obj,xNew] = obj.update();
                %% penalty factor continuation
                obj.m_solver = obj.m_solver.performPenaltyContinuation();
                
                %% check design change
                delta = abs(xNew(:) - obj.m_x(:));
                obj.m_change = max(delta(:));
                obj.m_x = xNew;
               
            end
            if (iter == obj.m_maxNumIters)
                obj.m_flag = 2; % optimization terminated after max iterations
                if (~obj.m_testMode),disp(['Optimization terminated after ' num2str(iter)  ' iterations']);end
            end
            if (~obj.m_testMode)
                obj= obj.filterDensity();
                obj.m_solver = obj.m_solver.setDesign(obj.m_xPhys);
                obj.plotPseudoDensity();
                obj.plotIsoSurface('Contour');
                obj.plotConvergence();
            end
        end
        %% SOLVE FEA AND ADJOINT
        function obj = solve(obj)
            obj.m_solver = obj.m_solver.setDesign(obj.m_xPhys);
            obj.m_solver = obj.m_solver.solve();
            obj.m_solver = obj.m_solver.postProcess();
        end
        %% CHECK TERMINATION
        function [obj,isConverged] = checkConvergence(obj)
            % Check for convergence
            nCovergenceIters = 10;
            isConverged = false;
            if (obj.m_iter > nCovergenceIters && ...
                    abs(obj.m_history.constraint.volFrac(obj.m_iter)-obj.m_targetVolFrac) < 0.02)
                VolChange = max(abs(obj.m_history.constraint.volFrac(obj.m_iter)-obj.m_history.constraint.volFrac(obj.m_iter-nCovergenceIters:obj.m_iter-1)));
                ObjChange = max(abs(obj.m_history.objective(obj.m_iter)-obj.m_history.objective(obj.m_iter-nCovergenceIters:obj.m_iter-1)))/obj.m_history.objective(1);
                DensityChange = max(obj.m_history.change(obj.m_iter-nCovergenceIters:obj.m_iter));
                if (VolChange < 0.02) && (ObjChange < 0.02) && (DensityChange < obj.m_changeTarget)
                    isConverged = true;
                end
            end
        end
        %% EVALUATE OVERRIDE
        function obj = evaluate(obj)
            % evaluate objective
            [obj.m_objective,obj.m_fx] = obj.m_objective.evaluate();
            obj.m_fx = mean(obj.m_fx); % take mean if needed for number of scenarios
            if (obj.m_iter==1)
                obj.m_fx0 = obj.m_fx;  
                obj.m_scale = obj.m_objective.m_scale;
            end

            % evaluate constraints
            obj.m_gx = zeros(obj.m_numConstraints,1);
            for g = 1 : obj.m_numConstraints
                [obj.m_constraints{g}, gx] = obj.m_constraints{g}.evaluate();
                obj.m_gx(g) = gx;
            end
        end
        %% GRADIENT OVERRIDE
        function obj = gradient(obj)
            % gradient objective
            obj.m_dfdx = zeros(obj.m_solver.m_ny,obj.m_solver.m_nx,obj.m_solver.m_numScenarios);
            [obj.m_objective,obj.m_dfdx] = obj.m_objective.gradient();
            obj.m_dfdx = mean(obj.m_dfdx,3);
            obj.m_dfdx = obj.m_dfdx(:);
            % gradient constraints
            obj.m_dgdx = [];
            for g = 1 : obj.m_numConstraints
                [obj.m_constraints{g}, dgdx] = obj.m_constraints{g}.gradient();
                dgdx = dgdx(:); % convert to vector
                dgdx = dgdx';
                obj.m_dgdx = [obj.m_dgdx;dgdx]; % append
            end
        end
        %% UPDATE DESIGN
        function [obj,xnew] = update(obj)
            if (strcmp(obj.m_update,'OC') == 1)
                [obj,xnew] = obj.OC();
            elseif (strcmp(obj.m_update,'MMA') == 1)
                [obj,xnew] = obj.MMA();
            elseif (strcmp(obj.m_update,'GCMMA') == 1)
                [obj,xnew] = obj.GCMMA();
            end
        end
        function [obj,xnew] = OC(obj)
            dfdx = reshape(obj.m_dfdx, obj.m_solver.m_ny, obj.m_solver.m_nx);
            [obj.m_ocOptimizer,xnew] = obj.m_ocOptimizer.update(obj.m_x,dfdx);
        end
        function [obj,xNew]=MMA(obj)
            xmma = obj.m_x(:);
            xmax = ones(size(xmma));
            xmin =  1e-9*xmax;
            fx = obj.m_fx/obj.m_scale;
            dfdx = obj.m_dfdx(:)/obj.m_scale;
            gx = obj.m_gx;
            dgdx = obj.m_dgdx;
            [obj.m_mmaOptimizer,xmma] = obj.m_mmaOptimizer.update( ...
                xmma,fx,dfdx,gx,dgdx,xmin,xmax);
            xNew = reshape(xmma,size(obj.m_x));
            xNew = xNew .* obj.m_solver.m_existingElems;
        end
        function [obj,xNew]=GCMMA(obj)
            xmma = obj.m_x(:);
            xmax = ones(size(xmma));
            xmin = 1e-9*ones(size(xmma));
            % SOLVE OUTER ITERATION
            fx = obj.m_fx/obj.m_scale;
            dfdx = obj.m_dfdx(:)/obj.m_scale;
            gx = obj.m_gx;
            dgdx = obj.m_dgdx;
            [obj.m_gcmmaOptimizer,xnew_mma] = obj.m_gcmmaOptimizer.outerUpdate(xmma,fx,dfdx,gx,dgdx,xmin,xmax);
            obj.m_x = reshape(xnew_mma,size(obj.m_x));
            % filter density
            obj= obj.filterDensity();
            obj.m_solver = obj.m_solver.setDesign(obj.m_xPhys);

            obj = obj.solve();
            obj = obj.evaluate();
            fxnew = obj.m_fx/obj.m_scale;
            gxnew = obj.m_gx;
            conserv = obj.m_gcmmaOptimizer.conCheck(fxnew,gxnew);
            % SOLVE INNERE ITERATIONS
            innerit=0;
            if conserv == 0
                while conserv == 0 && innerit <= 5
                    innerit = innerit+1;
                    [obj.m_gcmmaOptimizer,xnew_mma] = obj.m_gcmmaOptimizer.innerUpdate(xnew_mma, ...
                        fxnew,gxnew, xmma,fx, dfdx, gx, dgdx, xmin,xmax);
                    obj.m_x = reshape(xnew_mma,size(obj.m_x));
                    % filter density
                    obj= obj.filterDensity();

                    obj.m_solver = obj.m_solver.setDesign(obj.m_xPhys);
                    obj = obj.solve();
                    obj = obj.evaluate();
                    fxnew = obj.m_fx/obj.m_scale;
                    gxnew = obj.m_gx;
                    conserv = obj.m_gcmmaOptimizer.conCheck(fxnew,gxnew);
                end
            end
            % Update xolds
            obj.m_gcmmaOptimizer.m_xold2 = obj.m_gcmmaOptimizer.m_xold1;
            obj.m_gcmmaOptimizer.m_xold1 = xmma;
            obj.m_x = reshape(xmma,size(obj.m_x));
            xNew = reshape(xnew_mma,size(obj.m_x));
        end

        %% FILTERS
        function obj = filterDensity(obj)
            obj.m_stageDesign = cell(obj.m_numMfgConstraints+1,1);
            obj.m_stageDesign{1} = obj.m_x .* obj.m_solver.m_existingElems;  % x^(0)

            for k = 1:obj.m_numMfgConstraints
                obj.m_stageDesign{k+1} = obj.m_mfgConstraints{k}.filterDesign(obj.m_stageDesign{k});
            end
            obj.m_xPhys = obj.m_stageDesign{end}; % x^(M)
        end

        function obj = filterSensitivity(obj)
            numElems = obj.m_solver.m_nx*obj.m_solver.m_ny;
            num_fIds = round(length(obj.m_dfdx)/numElems);
            num_gIds = round(length(obj.m_dgdx)/numElems);

            % Objectives
            for fId = 1:num_fIds
                dfdx = obj.m_dfdx((fId-1)*numElems+1:fId*numElems,1);
                dfdx = reshape(dfdx, obj.m_solver.m_ny, obj.m_solver.m_nx);  % df/dx^(M)
                sens = dfdx;
                for k = obj.m_numMfgConstraints:-1:1
                    design_in = obj.m_stageDesign{k};  % x^(k-1)
                    sens = obj.m_mfgConstraints{k}.filterSensitivity(design_in, sens);
                end
                sens = sens .* obj.m_solver.m_existingElems;
                obj.m_dfdx((fId-1)*numElems+1:fId*numElems,1) = sens(:);
            end
            % Constraints
            for gId = 1:num_gIds
                % Read as column slice (df/dx^(M) for constraint gId)
                dgdx = obj.m_dgdx(gId,:);
                dgdx = dgdx(:);
                dgdx = reshape(dgdx, obj.m_solver.m_ny, obj.m_solver.m_nx);

                % Backward chain through manufacturing constraints
                sens = dgdx .*obj.m_solver.m_existingElems;
                for k = obj.m_numMfgConstraints:-1:1
                    x_in = obj.m_stageDesign{k};
                    sens = obj.m_mfgConstraints{k}.filterSensitivity(x_in, sens);
                end
                dgdx = sens(:);
                obj.m_dgdx(gId,:) = dgdx';
            end
        end
        %% OUTPUT AND VISUALIZATION
        function obj = plotConvergence(obj)
            plt = PlotId;
            figure(plt.convergence); set(gcf, 'Name', 'Convergence')
            plot(1:obj.m_iter,obj.m_history.objective(1:obj.m_iter)/obj.m_history.objective(1),'-b','LineWidth',2);
            xlabel('Iteration');
            if (strcmp(obj.m_objective,'compliance') == 1)
                ylabel('Normalized Compliance $C/C_0$');
            elseif (strcmp(obj.m_objective,'vonmises'))
                ylabel('Normalized Stress P-norm $\sigma^{PN}/\sigma^{PN}_0$');
            elseif (strcmp(obj.m_objective,'displacement') == 1)
                ylabel('Normalized Displacement $d/d_0$');
            elseif (strcmp(obj.m_objective,'dissipatedEnergy') == 1)
                ylabel('Normalized Dissipated Energy $\Xi/\Xi_0$');
            else
                ylabel('Normalized Objective $\varphi/\varphi_0$');
            end
            axis tight;
            ylim([0 inf])
            yyaxis right
            plot(1:obj.m_iter,obj.m_history.constraint.volFrac(1:obj.m_iter),'--r','LineWidth',2);
            ylabel('Volume Fraction ($V/V_0$)');
            ylim([0,1])
        end

        function plotPseudoDensity(obj)
            plt = PlotId;
            cm = ColorMaps;
            fig = figure(plt.design); clf(fig,'reset');
            set(gcf, 'Name', 'Pseudo-density');
            imagesc(1-flipud(obj.m_xPhys));
            colormap(cm.design);
            pbaspect(obj.m_solver.m_boxSizes);axis off; grid off;axis tight;
            pause(1e-4)
            drawnow
        end

        function plotIsoSurface(obj,method)
            plt = PlotId;
            cm = ColorMaps;
            fig = figure(plt.isosurface_contour);clf(fig,'reset');
            set(gcf, 'Name', 'Isosurface');
            if (strcmp(method,'Contour') == 1)
                colormap(cm.isosurface_contour)
                contourf(-obj.m_xPhys,[-0.5 -1]);
            else
                disp(['Method ' method 'for iso-surface generation is not implemented!']);
            end
            axis off
            pbaspect(obj.m_solver.m_boxSizes);axis off; grid off;axis tight;
            drawnow
        end
    end
end