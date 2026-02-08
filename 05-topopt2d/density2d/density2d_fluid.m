%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                              %
% This class implements a 2D density-based topology optimization using      %
% fluid flow problems. It inherits from the density2d class and implements  %
% the methods necessary to solve the fluid flow problem.                    %       
% % This code is largely based on the MATLAB code written by Joe            %
% Alexandersen                                                              %
% included in the 'utilities/thirdParty/topflow' directory. For more        %
% details, please refer to the original paper:                              %
% Alexandersen, Joe. "A detailed introduction to density-based              %
% topology optimisation of fluid flow problems with                         %
% implementation in MATLAB." Structural and Multidisciplinary               %
% Optimization 66.1 (2023): 12.                   %
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

classdef density2d_fluid < density2d
    properties
        % dissipatedEnergy minimization
        m_cx0; % initial dissipatedEnergy value
        m_cx; % dissipatedEnergy value
        m_dcdx; % sensitivity of dissipatedEnergy
    end

    methods
        function obj = density2d_fluid(solver,...
                objective,constraints,mfgConstraints, ...
                update, ...
                maxNumIters,exportGIF,testMode)

            % check inputs
            if ~isa(solver, 'fea2d_fluid')
                error('solver must be an instance of fea2d_fluid class');
            end
            % set default values
            if (nargin < 4)
                mfgConstraints.rmin = 1.5;
            end
            if (nargin < 5),update = 'MMA';end
            if (nargin < 6),maxNumIters = 250;end
            if (nargin < 7),exportGIF = false;end
            if (nargin < 8),testMode = false;end

            % construct
            obj = obj@density2d(solver,objective,constraints,mfgConstraints, ...
                update, ...
                maxNumIters,exportGIF,testMode); % call superclass
            
            obj.m_x = ones(size(obj.m_solver.m_existingElems));
            obj.m_x(obj.m_solver.m_activeDesignDomain==1) = obj.m_targetVolFrac;
            obj.m_xPhys = obj.m_x;
            obj.m_solver = obj.m_solver.setDesign(obj.m_x);

            obj.m_solver = obj.m_solver.setupContinuationScheme();
            % history
            obj.m_history.constraint.volFrac = zeros(obj.m_maxNumIters,1);
            obj.m_history.state.velocity = zeros(obj.m_maxNumIters,obj.m_solver.m_numScenarios);
            obj.m_history.state.pressure = zeros(obj.m_maxNumIters,obj.m_solver.m_numScenarios);
        end        

        function obj = set_qa(obj,qa)
            obj.m_solver = obj.m_solver.set_qa(qa);
        end
        %% FILTERS
        function obj = filterDensity(obj)
            
            idx = obj.m_solver.m_activeDesignDomain==1;
            obj.m_x(~idx) = 1;
            obj.m_stageDesign = cell(obj.m_numMfgConstraints+1,1);
            obj.m_stageDesign{1} = obj.m_x .* obj.m_solver.m_existingElems;  % x^(0)

            for k = 1:obj.m_numMfgConstraints
                obj.m_stageDesign{k+1} = obj.m_mfgConstraints{k}.filterDesign(obj.m_stageDesign{k});
            end
            obj.m_xPhys(idx) = obj.m_stageDesign{end}(idx); % x^(M)
        end

        function obj = filterSensitivity(obj)
            idx = obj.m_solver.m_activeDesignDomain==1;
            numElems = obj.m_solver.m_nx*obj.m_solver.m_ny;
            num_fIds = round(length(obj.m_dfdx)/numElems);
            num_gIds = round(length(obj.m_dgdx)/numElems);

            % Objectives
            for fId = 1:num_fIds
                dfdx = obj.m_dfdx((fId-1)*numElems+1:fId*numElems,1);
                dfdx = reshape(dfdx, obj.m_solver.m_ny, obj.m_solver.m_nx);  % df/dx^(M)

                sens = dfdx .* obj.m_solver.m_activeDesignDomain;
                for k = obj.m_numMfgConstraints:-1:1
                    design_in = obj.m_stageDesign{k};  % x^(k-1)
                    sens = obj.m_mfgConstraints{k}.filterSensitivity(design_in, sens);
                end
                sens = sens .* obj.m_solver.m_activeDesignDomain;
                obj.m_dfdx((fId-1)*numElems+1:fId*numElems,1) = sens(:);
            end
            % Constraints
            for gId = 1:num_gIds
                % Read as column slice (df/dx^(M) for constraint gId)
                dgdx = obj.m_dgdx(gId,:);
                dgdx = dgdx(:);
                dgdx = reshape(dgdx, obj.m_solver.m_ny, obj.m_solver.m_nx);

                % Backward chain through manufacturing constraints
                sens = dgdx .*obj.m_solver.m_activeDesignDomain;
                sens(~idx) = 0;
                for k = obj.m_numMfgConstraints:-1:1
                    x_in = obj.m_stageDesign{k};
                    sens = obj.m_mfgConstraints{k}.filterSensitivity(x_in, sens);
                end
                sens(~idx) = 0;
                dgdx = sens(:);
                obj.m_dgdx(gId,:) = dgdx';
            end
        end
        %% HISTORY
        function obj = saveHistory(obj)
            idx = obj.m_solver.m_activeDesignDomain==1;
            x_active = obj.m_xPhys(idx);
            numActiveDomain = sum(obj.m_solver.m_activeDesignDomain(:));
            obj.m_history.constraint.volFrac(obj.m_iter) = 1-sum(x_active(:))/numActiveDomain;
            obj.m_history.objective(obj.m_iter)= obj.m_fx;
            obj.m_history.change(obj.m_iter) = obj.m_change;
            for scenarioId = 1:obj.m_solver.m_numScenarios
                obj.m_history.state.velocity(obj.m_iter,scenarioId) = max(obj.m_solver.m_velocity.norm(:,:,scenarioId),[],'all');
                obj.m_history.state.pressure(obj.m_iter,scenarioId) = max(obj.m_solver.m_pressure(:,:,scenarioId),[],'all');
            end
        end
        function obj = setPseudoDensityInRectangle(obj,center,w,h,rho_in,rho_out)
            % setPseudoDensityInRectangle  Mark elements inside a rectangle with pseudo-density = 1
            %
            %   center : [cx, cy] center of rectangle (in global coordinates)
            %   w, h   : full width and height of rectangle
            %
            %   The rectangle is axis-aligned and elements are marked based on their
            %   center coordinates in obj.m_elemCoords. Only existing elements
            %   (m_existingElems ~= 0) are considered.

            cx = center(1);
            cy = center(2);

            xMinRect = cx - 0.5*w;
            xMaxRect = cx + 0.5*w;
            yMinRect = cy - 0.5*h;
            yMaxRect = cy + 0.5*h;

            % Element center coordinates as column vectors
            xe = obj.m_solver.m_elemCoords(1, :).';
            ye = obj.m_solver.m_elemCoords(2, :).';

            % Mask of existing elements in linear indexing
            maskExisting = obj.m_solver.m_existingElems(:) ~= 0;

            % Mask of elements whose centers are inside the rectangle
            insideRect = (xe >= xMinRect) & (xe <= xMaxRect) & ...
                (ye >= yMinRect) & (ye <= yMaxRect);

            % Combine masks and write into pseudo-density field
            idx = find(maskExisting & insideRect);
            obj.m_x = rho_out*ones(size(obj.m_solver.m_existingElems));
            obj.m_x(idx) = rho_in;

            obj.m_solver = obj.m_solver.setDesign(obj.m_x);
        end
        %% OUTPUT
        function obj = printResults(obj)

            disp(['Iteration: ' sprintf('%4i',obj.m_iter) ...
                ', Obj.: ' sprintf('%1.2e',obj.m_history.objective(obj.m_iter)) ...
                ', V/V0: ' sprintf('%1.2f',obj.m_history.constraint.volFrac(obj.m_iter))...
                ', f/f0: ' sprintf('%1.2f',obj.m_history.objective(obj.m_iter)/obj.m_history.objective(1)) ...
                ', Design Change: ' sprintf('%1.4f',obj.m_history.change(obj.m_iter))]);

            for scenarioId = 1:obj.m_solver.m_numScenarios
                disp(['Scenario: ' sprintf('%4i',scenarioId) ...
                    ', Velocity: ' sprintf('%1.2e',obj.m_history.state.velocity(obj.m_iter,scenarioId)) ...
                    ', Pressure: ' sprintf('%1.2e',obj.m_history.state.pressure(obj.m_iter,scenarioId)) ...
                    ])
            end
            fprintf('\n');
        end
        %% PLOT OVERRIDE
        function plotPseudoDensity(obj)
            plt = PlotId;
            cm = ColorMaps;
            fig = figure(plt.design); clf(fig,'reset');
            set(gcf, 'Name', 'Pseudo-density');
            imagesc(flipud(obj.m_xPhys));
            colormap(cm.design);
            pbaspect(obj.m_solver.m_boxSizes);axis off; grid off;axis tight;
            drawnow
            pause(1e-4)
        end

        function plotIsoSurface(obj,method)
            plt = PlotId;
            cm = ColorMaps;
            fig = figure(plt.isosurface_contour);clf(fig,'reset');
            set(gcf, 'Name', 'Isosurface');
            if (strcmp(method,'Contour') == 1)
                colormap(cm.isosurface_contour)
                contourf(obj.m_xPhys,[0 0.5]);
            else
                disp(['Method ' method 'for iso-surface generation is not implemented!']);
            end
            pbaspect(obj.m_solver.m_boxSizes);axis off;axis tight; grid off;
            drawnow
            pause(1e-4)
        end
    end
end