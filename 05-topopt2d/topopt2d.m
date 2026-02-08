%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description:                                                              %
% This class implements an abstract class for 2D topology optimization.     %
% It desfines the API for the optimization process.                         %
% It also provides common methods for retaining elements,                   %
% dilating the retained elements, and plotting convergence.                 %
%                                                                           %
% The class is intended to be inherited by specific topology optimization   %
% methods such as density, level-set, evolutionary, standard and            %
% modified Hamilton-Jacobi, and Pareto-tracing methods.                     %
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

classdef (Abstract) topopt2d < handle
    properties
        m_flag = 0; % flag for the optimization process, 0: not started, -1: failed, 1: converged, 2: max iterations reached, 3: user stopped
        m_solver; % fea solver
        m_x; % design variables (pseudo-densities)
        m_xPhys; % physical pseudo-densities
        m_objective; % objective type
        m_constraints; % list of constriants and upper bounds
        m_numConstraints; % number of constraints

        m_fx0; % initial objective function
        m_fx; % objective function
        m_gx; % constraints
        m_dfdx; % objective sensitivity
        m_dgdx; % constraint sensitivity

        m_mfgConstraints; % manufacturing constriants
        m_numMfgConstraints; % number of manufacturing constriants
        m_retainElems; % elements to be retained

        m_history; % optimization history
        m_iter; % current iteration
        m_maxNumIters; % maximum number of iterations
        m_change; % change in density w.r.t. previous step
        m_changeTarget = 0.005; % density change for termination

        m_retainNeighborIds; % for retain

        m_upsampledDesign; % upsampled for plot and DXF export (nyF x nxF), contains NaNs outside
        m_upsampled_xF;      % 1 x (nxF)
        m_upsampled_yF;      % 1 x (nyF)

        m_exportGIF = false; % export .gif file of the optimization process
        m_testMode = false; % test mode flag, if true, the optimization will not run and only the initial setup will be done
    end
    methods (Abstract)
        obj = optimize(obj)
        obj = solve(obj)
        obj = update(obj)
        obj = evaluate(obj)
        obj = gradient(obj)
    end
    methods
        function obj = topopt2d(solver,objective,constraints, mfgConstraints, ...
                maxNumIters,exportGIF, testMode)

            obj.m_solver = solver;
            obj.m_objective = objective;
            obj.m_constraints = constraints;
            obj.m_mfgConstraints = mfgConstraints;
            obj.m_numMfgConstraints = length(obj.m_mfgConstraints);

            obj.m_maxNumIters = maxNumIters;
            obj.m_change = 0;

            % retain neighborhood
            obj.m_retainElems = zeros(size(obj.m_solver.m_existingElems));
            obj.m_retainNeighborIds = [
                0  0;  % Center (current element)
                0 -1;  % Left
                -1  0;  % Below
                0  1;  % Right
                1  0;  % Above
                ];

            obj.m_numConstraints = length(obj.m_constraints); % number of constraints

            obj.m_exportGIF = exportGIF;
            obj.m_testMode = testMode;

            %% History
            % objective: Assuming 1 objective function (e.g., average of
            % compliance for all load scenarios)
            %
            % change: Fraction of change in design ( e.g., change in
            % volume)
            %
            % constraint: structure for all relevant constraints (e.g., volume)
            %
            % state: maximum values of the state variables and other
            % derived fields (temperature, deformation, von Mises stress)
            % at each load scenario

            obj.m_history = struct( ...
                'objective',zeros(obj.m_maxNumIters,1),...
                'change',zeros(obj.m_maxNumIters,1), ... %
                'constraint',struct(), ... % Constraints structure (e.g., volume)
                'state', struct() ... %
                );

        end

        function retainElems = dilateRetain(obj,retainElems)
            r = max(1,floor(obj.m_mfgConstraints.rmin));
            se = strel('disk', r);
            retainElems = imdilate(retainElems, se);
        end

        function obj = retainEdge(obj,seg)
            retainElems = zeros(obj.m_solver.m_ny,obj.m_solver.m_nx);
            nNode_y = obj.m_solver.m_ny+1;
            for segId = 1:numel(seg)
                segNodes = obj.m_solver.findNodesOnEdge(seg(segId));
                for k = 1:numel(segNodes)
                    node = segNodes(k);
                    yId = mod(node,nNode_y);
                    if yId == 0
                        yId = nNode_y; % Correct for MATLAB's 1-based indexing
                    end
                    xId = ceil((node-yId)/nNode_y);
                    for neighborId = 1:size(obj.m_retainNeighborIds,1)
                        j = yId + obj.m_retainNeighborIds(neighborId,1);
                        i = xId + obj.m_retainNeighborIds(neighborId,2);
                        neighbor = [j,i];
                        if (neighbor(1)>0 && neighbor(1)<= obj.m_solver.m_ny && ...
                                neighbor(2)>0 && neighbor(2)<= obj.m_solver.m_nx)
                            if (obj.m_solver.m_existingElems(neighbor(1),neighbor(2))==1)
                                retainElems(neighbor(1),neighbor(2)) = 1;
                            end
                        end
                    end
                end
            end
            % dilate
            retainElems = obj.dilateRetain(retainElems);
            obj.m_retainElems = obj.m_retainElems | (retainElems==1);
            obj.m_retainElems = obj.m_retainElems .* obj.m_solver.m_existingElems;
        end

        function obj = retainNeumannElements(obj)
            retainElems = zeros(obj.m_solver.m_ny,obj.m_solver.m_nx);
            nNode_y = obj.m_solver.m_ny + 1;
            for scenarioId = 1:obj.m_solver.m_numScenarios
                scenarioNodes = obj.m_solver.m_forcedNodes{scenarioId};
                for k = 1:numel(scenarioNodes)
                    node = scenarioNodes(k);
                    yId = mod(node,nNode_y);
                    if yId == 0
                        yId = nNode_y; % Correct for MATLAB's 1-based indexing
                    end
                    xId = ceil((node-yId)/nNode_y);
                    for neighborId = 1:size(obj.m_retainNeighborIds,1)
                        j = yId + obj.m_retainNeighborIds(neighborId,1);
                        i = xId + obj.m_retainNeighborIds(neighborId,2);
                        neighbor = [j,i];
                        if (neighbor(1)>0 && neighbor(1)<= obj.m_solver.m_ny && ...
                                neighbor(2)>0 && neighbor(2)<= obj.m_solver.m_nx)
                            if (obj.m_solver.m_existingElems(neighbor(1),neighbor(2))==1)
                                retainElems(neighbor(1),neighbor(2)) = 1;
                            end
                        end
                    end
                end
            end
            retainElems = obj.dilateRetain(retainElems);
            obj.m_retainElems = obj.m_retainElems | (retainElems==1);
            obj.m_retainElems = obj.m_retainElems .* obj.m_solver.m_existingElems;
        end


        function obj = upsampleDesignField(obj, field)
            % element-centered fields
            phiE  = field;                                 % (ny x nx)
            maskE = (obj.m_solver.m_existingElems == 1); % (ny x nx)

            % Do NOT NaN-clip with maskE here (causes jagged domain boundary)

            [ny,nx] = size(phiE);

            % 1) Smooth mask boundary via signed-distance field on elements
            % inside negative, outside positive
            din  = bwdist(~maskE);
            dout = bwdist(maskE);
            sdE  = dout - din;  % ~0 at mask boundary

            % 2) Convert element-centered -> nodal by averaging (phi + sd)
            phiN = nan(ny+1,nx+1);
            sdN  = nan(ny+1,nx+1);

            for j = 1:ny+1
                js = max(1,j-1):min(ny,j);
                for i = 1:nx+1
                    is = max(1,i-1):min(nx,i);

                    vphi = phiE(js,is); vphi = vphi(~isnan(vphi));
                    if ~isempty(vphi)
                        phiN(j,i) = mean(vphi);
                    end

                    vsd = sdE(js,is);
                    sdN(j,i) = mean(vsd(:));
                end
            end

            % 3) Coordinates in physical units (nodal grid)
            Xmin = obj.m_solver.m_boundingBox(1,1);
            Xmax = obj.m_solver.m_boundingBox(1,2);
            Ymin = obj.m_solver.m_boundingBox(2,1);
            Ymax = obj.m_solver.m_boundingBox(2,2);

            xN = linspace(Xmin, Xmax, nx+1);
            yN = linspace(Ymin, Ymax, ny+1);
            [XN,YN] = meshgrid(xN,yN);

            % 4) Upsample with LINEAR interpolation (phi + mask SDF)
            ref = 6; % increase if needed
            xF  = linspace(xN(1), xN(end), ref*nx + 1);
            yF  = linspace(yN(1), yN(end), ref*ny + 1);
            [XF,YF] = meshgrid(xF,yF);

            phiF = interp2(XN,YN,phiN,XF,YF,'linear');
            sdF  = interp2(XN,YN,sdN, XF,YF,'linear');

            % 5) Clip by SMOOTH mask boundary (prevents jagged domain edge)
            phiF(sdF > 0) = NaN;

            obj.m_upsampledDesign = phiF;    % (nyF x nxF), contains NaNs outside
            obj.m_upsampled_xF    = xF;      % 1 x (nxF)
            obj.m_upsampled_yF    = yF;      % 1 x (nyF)
        end

        %% OUTPUT AND VISUALIZATION
        function obj = plotConvergence(obj)
            % PLOT CONVERGENCE
            plt = PlotId;
            figure(plt.convergence); set(gcf, 'Name', 'convergence')
            plot(1:obj.m_iter,obj.m_history.objective(1:obj.m_iter));

            xlabel('Iteration');
            if (strcmp(obj.m_objective,'compliance') == 1),ylabel('Relative Compliance')
            else, ylabel('Objective'); end
        end

        function obj = exportDXF(obj, filename)
            % PLOT AND EXPORT DXF LINES
            if (isempty(obj.m_upsampledDesign))
                obj = obj.upsampleDesignField(obj.m_xPhys);
            end

            phiF = -obj.m_upsampledDesign;
            xF   = obj.m_upsampled_xF;
            yF   = obj.m_upsampled_yF;
            export_dxf_from_Levelset(phiF, xF, yF, -0.5, filename, 5);
        end

        function obj = exportSTL(obj, example_name,thickness)
            if nargin < 3; thickness = 1; end
            dxf_filename = [example_name '.dxf'];
            obj.exportDXF(example_name);
            dxf2stl(dxf_filename,thickness);
        end
    end
end