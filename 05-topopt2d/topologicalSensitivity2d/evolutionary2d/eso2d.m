%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                              %
% This class implements a 2D evolutionary-based topology optimization using %
% evolutionary optimization. It inherits from the topopt2d class and        %
% implements the methods necessary to solve the evolutionary structural     %
% optimization (ESO) problems.                                              %
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

classdef (Abstract) eso2d  < topopt2d
    properties
        m_targetVolFrac; % targt volume fraction
        m_volDecrement; % step size for material removal
        m_vx; % volume
        m_tau; % levelset threshold
    end
    methods (Abstract)
        obj = evaluate(obj)
        obj = gradient(obj)
        obj = saveHistory(obj)
    end
    methods
        function obj = eso2d(solver,objective,constraints, mfgConstraints, ...
                volDecrement,exportGIF,testMode)

            % construct
            maxNumIters = 300;
            obj = obj@topopt2d(solver,objective,constraints,mfgConstraints, ...
                maxNumIters,exportGIF,testMode); % call superclass

            for consId = 1 : obj.m_numConstraints
                if (isa(obj.m_constraints{consId},'volume'))
                    obj.m_targetVolFrac = obj.m_constraints{consId}.m_upperBound;
                    break;
                end
            end

            obj.m_volDecrement = volDecrement;

            obj.m_x = obj.m_solver.m_existingElems;
            obj.m_solver = obj.m_solver.setDesign(obj.m_x);
            obj.m_solver = obj.m_solver.setPenaltyFactor(1);
        end

        function obj = optimize(obj)
            obj.m_flag = -1; % optimization flag, intialize to failed
            obj.m_iter = 1; % current iteration
            
            obj = obj.evaluateVolume();
            volFrac = obj.m_vx; % current volume fraction

            obj.m_tau = 0;
            obj = obj.filterDensity();
            obj = obj.solve();
            obj = obj.evaluate();
            obj = obj.gradient();
            obj= obj.filterSensitivity();
            obj = obj.saveHistory();
            if (~obj.m_testMode)
                obj = obj.printResults();
                obj.plotIsoSurface();
            end
            while obj.m_vx > obj.m_targetVolFrac
                obj = obj.evaluateVolume();
                volDec = min(obj.m_volDecrement,obj.m_vx - obj.m_targetVolFrac);
                if (obj.m_exportGIF), export_gifs();end
                if (volDec < 0.005),obj.m_flag = 1; break;end
                if (obj.m_fx/obj.m_fx0 > 10),break;end % divergence
                volFrac = volFrac - volDec;
                obj = obj.update(volFrac); % remove elements below threshold
                obj = obj.filterDensity();
                obj = obj.solve();
                obj = obj.evaluate();
                obj = obj.gradient();
                obj= obj.filterSensitivity();
                %% history and visualize
                obj.m_iter = obj.m_iter + 1;
                obj.m_change = volDec;
                obj = obj.saveHistory();
                if (~obj.m_testMode)
                    obj = obj.printResults();
                    obj.plotIsoSurface();
                end
            end
        end
        %% SOLVE FEA AND ADJOINT
        function obj = solve(obj)
            obj.m_solver = obj.m_solver.solve();
            obj.m_solver = obj.m_solver.postProcess();
        end
        %% UPDATE DESIGN
        function obj = update(obj,volFrac)
            % Find the level-set value such that the contour has given vol fraction
            obj.m_tau = obj.findContourValueWithVolumeFraction(volFrac);
            index = find(obj.m_dfdx < obj.m_tau); % eliminate all elements less than this value
            obj.m_x = obj.m_solver.m_existingElems; % start with the full domain
            obj.m_x(ind2sub(size(obj.m_dfdx),index)) = 0; % remove elements
            obj.m_solver = obj.m_solver.setDesign(obj.m_x);
        end

        function  value = findContourValueWithVolumeFraction(obj,volfrac)
            nElemsToRemove = (obj.m_solver.m_nx*obj.m_solver.m_ny - obj.m_solver.m_numExistingElems) + ...
                round(obj.m_solver.m_numExistingElems*(1-volfrac));
            sortedField = sort(obj.m_dfdx(:));
            value = sortedField(nElemsToRemove);
        end

        function obj = evaluateVolume(obj)
            obj.m_vx = sum(obj.m_x(:))/obj.m_solver.m_numExistingElems;
        end

        %% FILTERS
        function obj = filterDensity(obj)
            % filter density
            for mfgConsId = 1:obj.m_numMfgConstraints
                obj.m_x = obj.m_mfgConstraints{mfgConsId}.filterDesign(obj.m_x);
            end
        end

        function obj = filterSensitivity(obj)
            for mfgConsId = 1:obj.m_numMfgConstraints
                % Filter sensitivity
                obj.m_dfdx = obj.m_mfgConstraints{mfgConsId}.filterSensitivity(obj.m_x,obj.m_dfdx);
                % Apply existing elements
                obj.m_dfdx = obj.m_dfdx .* obj.m_solver.m_existingElems;
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
            figure(plt.pareto_front); set(gcf, 'Name', 'Evolution')
            plot(obj.m_history.constraint.volFrac(1:obj.m_iter), ...
                obj.m_history.objective(1:obj.m_iter), ...
                '-ko', 'LineWidth',2,'MarkerFaceColor','r');
            xlabel('Volume Fraction');
            if (strcmp(obj.m_objective,'compliance') == 1)
                ylabel('Normalized Compliance $C/C_0$');
            else
                ylabel('Normalized Objective $\varphi/\varphi_0$');
            end
        end
        function plotIsoSurface(obj,method)
            plt = PlotId;
            cm = ColorMaps;
            lsf = obj.m_dfdx;
            if nargin<2, method ='Contour';end
            if (strcmp(method,'Contour') == 1)
                % generate smooth contour from LSF
                obj = obj.upsampleDesignField(lsf);
                % plot
                fig = figure(plt.isosurface_contour); clf(fig,'reset');
                set(gcf, 'Name', 'Isosurface');
                colormap(cm.isosurface_contour)
                ax = axes(fig); hold(ax,'on');
                axis(ax,'equal','tight','off');
                set(ax,'YDir','normal');
                phiF = obj.m_upsampledDesign;    % (nyF x nxF), contains NaNs outside
                XF = obj.m_upsampled_xF;      % 1 x (nxF)
                YF = obj.m_upsampled_yF;      % 1 x (nyF)
                contourf(ax, XF, YF, -phiF, [-max(lsf(:)) -obj.m_tau ], 'LineColor','none');
                hold(ax,'off');
                axis off;
            elseif (strcmp(method,'LS') == 1)
                fig = figure(plt.isosurface_lsf);clf(fig,'reset');
                set(gcf, 'Name', 'Level-set Field'); axis off
                [X,Y]=meshgrid(1:obj.m_solver.m_nx,1:obj.m_solver.m_ny);
                maxVal = max(lsf(:));
                surf(X,Y,1+lsf/maxVal); hold on
                surf(X,Y,1+obj.m_tau/maxVal*ones(size(lsf)),'facealpha',0.25);hold on
                set(gca,'ZTickLabel',[]);
                contourf(X,Y,lsf/maxVal,[obj.m_tau/maxVal 1]);hold on
                hold off;
                view([-0.6 -0.9 0.7]);
                axis on
            else
                disp(['Method ' method 'for iso-surface generation is not implemented!']);
            end
            pbaspect(obj.m_solver.m_boxSizes);axis tight; grid off;
            drawnow;
            pause(0.001);
        end


        function obj = exportDXF(obj, filename,minPts)
            % PLOT AND EXPORT DXF LINES
            if (isempty(obj.m_upsampledDesign))
                lsf = obj.m_lsf(2:end-1,2:end-1);
                obj = obj.upsampleDesignField(lsf);
            end

            phiF = obj.m_upsampledDesign;
            xF   = obj.m_upsampled_xF;
            yF   = obj.m_upsampled_yF;
            export_dxf_from_Levelset(-phiF, xF, yF, -obj.m_tau, filename, minPts);
        end

        function obj = exportSTL(obj, example_name,thickness,minPts)
            if nargin < 3; thickness = 1; end
            if nargin < 4; minPts = 5; end
            dxf_filename = [example_name '.dxf'];
            obj.exportDXF(dxf_filename,minPts);
            dxf2stl(dxf_filename,thickness);
        end

    end

    methods (Static)
        function phiF = removeHangingRegions(phiF, tau)
            % Keep only the largest connected solid component (phiF >= tau).
            % Disconnected "hanging" solids are set to NaN so contourf ignores them.
            %
            % Needs Image Processing Toolbox (bwconncomp).

            solid = isfinite(phiF) & (phiF >= tau);   % solid mask (ignore NaNs)
            if ~any(solid(:)), return; end

            CC = bwconncomp(solid, 8);                % 8-connectivity is typical for pixels
            if CC.NumObjects <= 1, return; end

            % find largest component
            sz = cellfun(@numel, CC.PixelIdxList);
            [~, k] = max(sz);

            keep = false(size(phiF));
            keep(CC.PixelIdxList{k}) = true;

            % remove hanging solids
            remove = solid & ~keep;
            phiF(remove) = NaN;
        end
    end

end

