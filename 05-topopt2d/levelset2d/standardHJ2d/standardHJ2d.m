%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                              %
% This class implements a 2D level-set based topology optimization using    %
% the standard Hamilton-Jacobi (HJ) method. It inherits from the topopt2d    %
% class.                                                                    %
%                                                                           %
% The level-set approach and code is based on Challis, Vivien J.            %
% "A discrete level-set topology optimization code written in Matlab."      %
% Structural and multidisciplinary optimization 41 (2010): 453-464.         %
%                                                                           %
% For higher order solution, the LSF re-initialization to SDF is achieved   %
% using the "ToolboxLS" developed by:                                       %
% Ian M. Mitchell (mitchell@cs.ubc.ca), 2004                                %
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

classdef standardHJ2d < topopt2d
    properties
        m_targetVolFrac; % targt volume fraction
        m_numReinit; % Re-initializaion frequency
        m_stepLength;
        m_lsf; % Level-set function
        m_sdf; % Signed distance function

        m_velocity; % boundary velocity
        m_dSdx; % shape sensitivity of compliance function

        m_alpha;
        m_lagLower;
        m_lagUpper;

        m_nHolesX; % number of initial holes along x direction
        m_nHolesY; % number of initial holes along y direction
        m_r0; % radius of initial holes 0 < r0 < 1
    end
    methods (Abstract)
        obj = evaluate(obj)
        obj = gradient(obj)
        obj = saveHistory(obj)
    end
    methods
        function obj = standardHJ2d(solver,objective,constraints, mfgConstraints, ...
                nHolesX,nHolesY,r0, ...
                numReinit,stepLength,maxNumIters,exportGIF,testMode)

            % construct
            obj = obj@topopt2d(solver,objective,constraints,mfgConstraints, ...
                maxNumIters,exportGIF,testMode); % call superclass

            % target volume fraction
            for consId = 1 : obj.m_numConstraints
                if (isa(obj.m_constraints{consId},'volume'))
                    obj.m_targetVolFrac = obj.m_constraints{consId}.m_upperBound;
                    break;
                end
            end

            % initialize
            obj.m_x = obj.m_solver.m_existingElems;
            obj.m_solver = obj.m_solver.setDesign(obj.m_x);
            obj.m_solver = obj.m_solver.setPenaltyFactor(1);
            % Hole initializtion
            obj.m_nHolesX = nHolesX;
            obj.m_nHolesY = nHolesY;
            obj.m_r0 = r0;
            % LS parameters
            obj.m_maxNumIters = maxNumIters;
            obj.m_numReinit = numReinit;
            obj.m_stepLength = stepLength;
        end
        %% OPTIMIZE
        function obj = optimize(obj)
            obj.m_flag = -1; % optimization flag, intialize to failed
            % obj = obj.retainNeumannElements();
            if (obj.m_nHolesX*obj.m_nHolesY > 0)
                obj = obj.initializeHoles(obj.m_nHolesX,...
                    obj.m_nHolesY,obj.m_r0);
            end

            obj.m_solver = obj.m_solver.setDesign(obj.m_x);

            obj = obj.reinit();
            for iter = 1:obj.m_maxNumIters
                obj.m_iter = iter;
                %% solve
                obj = obj.solve();
                obj = obj.evaluate();
                obj = obj.gradient();
                %% history and display
                obj = obj.saveHistory();
                if (~obj.m_testMode)
                    obj = obj.printResults();
                    obj.plotIsoSurface('Contour');
                    if (obj.m_exportGIF), export_gifs();end
                end
                [obj,isConverged] = obj.checkConvergence();
                if (isConverged),obj.m_flag = 1; obj = obj.reinit();break;end
                %% update
                obj = obj.update();
                obj.m_solver = obj.m_solver.setDesign(obj.m_x);
                %% Reinitialize level-set function
                if ~mod(iter,obj.m_numReinit)
                    obj = obj.reinit();
                end
            end
            if (iter == obj.m_maxNumIters)
                obj.m_flag = 2; % optimization terminated after max iterations
                disp(['Optimization terminated after ' num2str(iter)  ' iterations']);
            end
            if (~obj.m_testMode)
                obj.plotPseudoDensity();
                obj.plotIsoSurface('LS');
                obj = obj.plotConvergence();
            end
        end
        function [obj,isConverged] = checkConvergence(obj)
            % Check for convergence
            nCovergenceIters = 5;
            isConverged = false;
            CC = bwconncomp(obj.m_x); % number of connected components

            if (CC.NumObjects == 1 && obj.m_iter > nCovergenceIters && ...
                    abs(obj.m_history.constraint.volFrac(obj.m_iter)-obj.m_targetVolFrac) < 0.005)
                VolChange = max(abs(obj.m_history.constraint.volFrac(obj.m_iter)-obj.m_history.constraint.volFrac(obj.m_iter-nCovergenceIters:obj.m_iter-1)));
                ObjChange = max(abs(obj.m_history.objective(obj.m_iter)-obj.m_history.objective(obj.m_iter-nCovergenceIters:obj.m_iter-1)))/obj.m_history.objective(1);
                if (VolChange < 0.002) && (ObjChange < 0.002)
                    isConverged = true;
                end
            end
        end
        %% REINITIALIZATION OF LEVEL-SET FUNCTION
        function obj = reinit(obj)
            if (isempty(obj.m_lsf))
                strucFull = zeros(size(obj.m_x)+2);
                strucFull(2:end-1,2:end-1) = obj.m_x;
                % Use "bwdist" (Image Processing Toolbox)
                obj.m_lsf = (~strucFull).*(bwdist(strucFull)-0.5) - ...
                    strucFull.*(bwdist(strucFull-1)-0.5);
            end
            % Use "ToolboxLS", for higher order solution
            % 2004, Ian M. Mitchell (mitchell@cs.ubc.ca)
            obj.m_sdf= reinit2D(obj.m_lsf,obj.m_solver.m_hx,obj.m_solver.m_boundingBox,'low');
            obj.m_lsf = obj.m_sdf;
        end
        %% SOLVE FEA AND ADJOINT
        function obj = solve(obj)
            obj.m_solver = obj.m_solver.solve();
            obj.m_solver = obj.m_solver.postProcess();
        end
        %% UPDATE DESIGN
        function obj = updateLagrangeMultipliers(obj)
            % Initialize Lagrangian multipliers and scaling parameters
            if obj.m_iter == 1
                % At the first iteration, initialize the lower and upper bounds for
                % Lagrange multipliers and the alpha scaling factor.
                obj.m_lagLower = -0.01;    % Initial lower bound for Lagrangian multiplier
                obj.m_lagUpper = 1000;    % Initial upper bound for Lagrangian multiplier
                obj.m_alpha = 0.9;        % Scaling factor for Lagrangian updates
            else
                % Update Lagrangian multipliers based on the current volume constraint
                % and the target volume fraction
                obj.m_lagLower = obj.m_lagLower - 1 / obj.m_lagUpper * ...
                    obj.m_gx;
                % Ensure the upper multiplier doesn't grow too small; adjust using alpha
                obj.m_lagUpper = max(obj.m_alpha * obj.m_lagUpper, 1);
            end
        end

        function obj = update(obj)
            ve = 1 / obj.m_solver.m_numExistingElems;
            % update Lagrange multipliers for augmented Lagrangian
            obj = obj.updateLagrangeMultipliers();

            % combine sensitivities based on augmented Lagrangian
            obj.m_dgdx = obj.m_dgdx * 1 / obj.m_lagUpper;

            obj.m_dSdx = obj.m_dfdx - ve*obj.m_lagLower + ...
                ve * 1/obj.m_lagUpper*obj.m_gx;

            % Smooth/filter the sensitivities to ensure stability and prevent oscillations
            obj = obj.filterSensitivity();

            % Extend sensitivities using a zero border to avoid boundary issues
            obj.m_velocity = zeros(size(obj.m_dSdx) + 2); % Create a padded array
            obj.m_velocity(2:end-1, 2:end-1) = -obj.m_dSdx; % Assign negative sensitivity values to the central region

            % Choose the time step based on the CFL condition to ensure numerical stability
            dt = 1 / (obj.m_stepLength * max(abs(obj.m_dSdx(:))));

            % Evolve the level set function over a total time determined by the CFL condition
            for i = 1:obj.m_stepLength
                % Compute forward and backward differences along the x-direction
                dpx = circshift(obj.m_lsf, [0, -1]) - obj.m_lsf; % Forward difference in x
                dmx = obj.m_lsf - circshift(obj.m_lsf, [0, 1]);  % Backward difference in x

                % Compute forward and backward differences along the y-direction
                dpy = circshift(obj.m_lsf, [-1, 0]) - obj.m_lsf; % Forward difference in y
                dmy = obj.m_lsf - circshift(obj.m_lsf, [1, 0]);  % Backward difference in y

                % Update the level set function using an upwind scheme
                % This ensures stable evolution of the level set based on velocity.
                obj.m_lsf = obj.m_lsf - dt * min(obj.m_velocity, 0) .* ...
                    sqrt(min(dmx, 0).^2 + max(dpx, 0).^2 + min(dmy, 0).^2 + max(dpy, 0).^2) ...
                    - dt * max(obj.m_velocity, 0) .* ...
                    sqrt(max(dmx, 0).^2 + min(dpx, 0).^2 + max(dmy, 0).^2 + min(dpy, 0).^2);
            end
            % Save the current design
            xOld = obj.m_x;

            % Extract the new structure from the updated level set function
            lsf = obj.m_lsf(2:end-1, 2:end-1);
            obj.m_x = (lsf < 0);

            % Filter the density
            obj = obj.filterDensity();

            % Evaluate the change in design variables
            delta_x = abs(xOld(obj.m_solver.m_existingElems == 1) - obj.m_x(obj.m_solver.m_existingElems == 1));
            obj.m_change = sum(delta_x(:)) / obj.m_solver.m_numExistingElems; % Normalized change
        end

        function obj = initializeHoles(obj,hx,hy,r)
            % Generate initial holes using the function:
            % Z = cos(X)*cos(Y)
            pltId = PlotId;
            [X,Y] = meshgrid(1:obj.m_solver.m_nx,1:obj.m_solver.m_ny);
            Z = cos(2*X*hx*pi/obj.m_solver.m_nx).*cos(2*Y*hy*pi/obj.m_solver.m_ny)+(1-r);
            Z = Z.*obj.m_solver.m_existingElems;
            obj.m_x = zeros(obj.m_solver.m_ny,obj.m_solver.m_nx);
            obj.m_x(Z > 0) = 1;

            % Filter the density to ensure retain elements are not removed
            obj = obj.filterDensity();

            if(~obj.m_testMode)
                % Plot Z function
                Z(obj.m_solver.m_existingElems==0) = nan;
                fig = figure(pltId.initial_holes); clf(fig,'reset');colormap(gray)
                contourf(-Z,[-max(Z(:)) 0 ]); pbaspect(obj.m_solver.m_boxSizes);
                set(gcf, 'Name', 'Initial Holes');
            end
        end

        %% FILTERS
        function obj = filterDensity(obj)
            % filter density
            for mfgConsId = 1:obj.m_numMfgConstraints
                obj.m_x = obj.m_mfgConstraints{mfgConsId}.filterDesign(obj.m_x);
            end
        end

        function obj = filterSensitivity(obj)
            % filter sensitivity of objective and constraints
            % the inputs are assumed to be flattened colummn vectors
            % each number of the vector corresponds to an element of the mesh
            % each nx by ny matrix corresponds to a scenario that is reshpaed to a matrix
            % of size nx by ny and then filtered and flattened again
            % The retain filter is applied to the each sensitivity field here
            % as well as the ensuring that we only have values for existing elements

            for mfgConsId = 1:obj.m_numMfgConstraints
                % Augment Shape Sensitivity
                obj.m_dSdx = obj.m_mfgConstraints{mfgConsId}.filterSensitivity(obj.m_x,obj.m_dSdx);
                % Apply existing elements
                obj.m_dSdx = obj.m_dSdx .* obj.m_solver.m_existingElems;
            end
        end

        %% OUTPUT
        function plotPseudoDensity(obj)
            plt = PlotId;
            cm = ColorMaps;
            F = obj.m_x;
            fig = figure(plt.design); clf(fig,'reset');
            set(gcf, 'Name', 'Pseudo-density');
            imagesc(1-flipud(F));
            colormap(cm.design);
            pbaspect(obj.m_solver.m_boxSizes);axis off; grid off;axis tight;
            pause(1e-4)
        end
        function obj = plotConvergence(obj)
            plt = PlotId;
            figure(plt.convergence); set(gcf, 'Name', 'Convergence')
            plot(1:obj.m_iter,obj.m_history.objective(1:obj.m_iter)/obj.m_history.objective(1),'-b','LineWidth',2);
            xlabel('Iteration');
            if (strcmp(obj.m_objective,'compliance') == 1)
                ylabel('Normalized Compliance $C/C_0$');
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

        function plotSensitivity(obj)
            plt = PlotId;
            cm = ColorMaps;
            dsdx = obj.m_dSdx;
            fig = figure(plt.ls_dsdx);clf(fig,'reset');
            set(gcf, 'Name', 'Sensitivity');
            colormap(cm.ls_dsdx)
            dsdx(obj.m_solver.m_existingElems == 0) = nan;
            [X,Y]=meshgrid(1:obj.m_solver.m_nx,1:obj.m_solver.m_ny);
            X = X*obj.m_solver.m_hx;
            Y = Y*obj.m_solver.m_hy;
            surf(X,Y,dsdx);
            axis off;
            view(2);
            pbaspect(obj.m_solver.m_boxSizes);axis tight; grid off;
        end

        function plotVelocity(obj)
            plt = PlotId;
            cm = ColorMaps;
            velocity = obj.m_velocity(2:end-1, 2:end-1);
            fig = figure(plt.ls_velcity); clf(fig, 'reset');
            set(gcf, 'Name', 'Velocity');
            colormap(cm.ls_velosity)

            [X, Y] = meshgrid(1:obj.m_solver.m_nx, 1:obj.m_solver.m_ny);
            X = X * obj.m_solver.m_hx;
            Y = Y * obj.m_solver.m_hy;

            % Compute the gradient of the velocity field
            [vx, vy] = gradient(velocity, obj.m_solver.m_hx, obj.m_solver.m_hy);

            % Normalize the gradients
            gradMagnitude = sqrt(vx.^2 + vy.^2);
            vx = vx ./ (gradMagnitude + eps);
            vy = vy ./ (gradMagnitude + eps);

            % Overlay arrows showing the velocity gradient
            quiver(X, Y, vx, vy, 1, 'r'); % Adjust scale factor as needed
            hold on;

            velocity(obj.m_solver.m_existingElems == 0) = nan;

            % Plot the velocity field
            surf(X, Y, velocity);
            axis off;
            view(2);
            pbaspect(obj.m_solver.m_boxSizes); axis tight; grid off;
        end

        function plotLSF(obj)
            plt = PlotId;
            cm = ColorMaps;
            lsf = obj.m_lsf(2:end-1,2:end-1);
            fig = figure(plt.lsf);clf(fig,'reset');
            set(gcf, 'Name', 'LSF');
            colormap(cm.lsf)
            lsf(obj.m_solver.m_existingElems == 0) = nan;
            [X,Y]=meshgrid(1:obj.m_solver.m_nx,1:obj.m_solver.m_ny);
            X = X*obj.m_solver.m_hx;
            Y = Y*obj.m_solver.m_hy;
            surf(X,Y,lsf);
            axis off;
            view(2);
            pbaspect(obj.m_solver.m_boxSizes);axis tight; grid off;
        end

        function plotIsoSurface(obj,method)
            plt = PlotId;
            cm  = ColorMaps;
            lsf = obj.m_lsf(2:end-1,2:end-1);

            if (strcmp(method,'Contour') == 1)
                % generate smooth contour from LSF
                obj = obj.upsampleDesignField(lsf);
                % plot
                fig = figure(plt.isosurface_contour); clf(fig,'reset');
                set(fig,'Name','Isosurface');
                colormap(cm.isosurface_contour);
                ax = axes(fig); hold(ax,'on');
                axis(ax,'equal','tight','off');
                set(ax,'YDir','normal');

                phiF = obj.m_upsampledDesign;    % (nyF x nxF), contains NaNs outside
                XF = obj.m_upsampled_xF;      % 1 x (nxF)
                YF = obj.m_upsampled_yF;      % 1 x (nyF)
                contourf(ax, XF, YF, phiF, [-Inf 0], 'LineColor','none');
                hold(ax,'off');

            elseif (strcmp(method,'LS') == 1)
                fig = figure(plt.isosurface_lsf);clf(fig,'reset');
                set(gcf, 'Name', 'Level-set Field'); axis off
                [X,Y]=meshgrid(1:obj.m_solver.m_nx,1:obj.m_solver.m_ny);
                X = X*obj.m_solver.m_hx;
                Y = Y*obj.m_solver.m_hy;
                maxVal = max(abs(lsf(:)));
                lsftmp = lsf;lsftmp(obj.m_solver.m_existingElems == 0) = nan;
                surf(X,Y,lsftmp/maxVal+1); hold on
                contourf(X,Y,-lsftmp/maxVal,[0 0]);
                surf(X,Y,ones(size(lsf)),'facealpha',0.25);hold on
                set(gca,'ZTickLabel',[]);

                hold off;
                view([-0.6 -0.9 0.7]);
                axis on;
            else
                disp(['Method ' method 'for iso-surface generation is not implemented!']);
            end
            pbaspect(obj.m_solver.m_boxSizes);axis tight; grid off;
            drawnow
            pause(0.001);
        end

        function obj = exportDXF(obj, filename)
            % PLOT AND EXPORT DXF LINES
            if (isempty(obj.m_upsampledDesign))
                lsf = obj.m_lsf(2:end-1,2:end-1);
                obj = obj.upsampleDesignField(lsf);
            end

            phiF = obj.m_upsampledDesign;
            xF   = obj.m_upsampled_xF;
            yF   = obj.m_upsampled_yF;
            export_dxf_from_Levelset(phiF, xF, yF, 0, filename, 5);
        end

        function obj = exportSTL(obj, example_name,thickness)
            if nargin < 3; thickness = 1; end
            dxf_filename = [example_name '.dxf'];
            obj.exportDXF(dxf_filename);
            dxf2stl(dxf_filename,thickness);
        end

    end
end

