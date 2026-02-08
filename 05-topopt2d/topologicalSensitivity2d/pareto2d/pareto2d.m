%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Description:                                                            %
% This class implements a 2D pareto-tracing topology optimization using     %
% topological sensitivity fields and fixed-point iteration at intermediate  %
% volume fractions. In other words, it traces the Pareto front of objective %
% and volume fraction and the intermediate designs are also locally optimal.%
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

classdef (Abstract) pareto2d < eso2d
    properties

        m_paretoAggressiveness; % higher (upto 1), closer to Pareto curve, but may terminate early

    end
    methods (Abstract)
        obj = evaluate(obj)
        obj = gradient(obj)
        obj = saveHistory(obj)
    end
    methods
        function obj = pareto2d(solver,objective,constraints, mfgConstraints, ...
                volDecrement,paretoAggressiveness,exportGIF,testMode)

            % construct
            obj = obj@eso2d(solver,objective,constraints, mfgConstraints, ...
                volDecrement,exportGIF,testMode) % call superclass

            obj.m_paretoAggressiveness = paretoAggressiveness;
        end

        %% UPDATE DESIGN
        function obj = update(obj,volFrac)
            iter = 0;
            isParetoOptimal = 0;
            while (iter < 20) % to avoid cycles typically a few iterations is sufficient
                if ((iter > 0)&&(isParetoOptimal)) % done with current vol
                    break
                end
                % Find the level-set value such that the contour has given vol fraction
                obj.m_tau = obj.findContourValueWithVolumeFraction(volFrac);
                index = find(obj.m_dfdx < obj.m_tau); % eliminate all elements less than this value
                obj.m_x = obj.m_solver.m_existingElems; % start with the full domain
                obj.m_x(ind2sub(size(obj.m_dfdx),index)) = 0; % remove elements
                obj.m_solver = obj.m_solver.setDesign(obj.m_x);
                obj = obj.filterDensity();
                obj = obj.solve();
                obj = obj.evaluate();
                obj = obj.gradient();
                obj= obj.filterSensitivity();
                isParetoOptimal = obj.analyzeTopology();
                iter= iter+1;
            end
        end

        function isParetoOptimal = analyzeTopology(obj)
            T_InMin = min(obj.m_dfdx(obj.m_x==1)); % Min of topological field inside the domain
            T_OutMax = max(obj.m_dfdx(obj.m_x==0 & obj.m_solver.m_existingElems == 1)); % Max of topological field outside the domain
            if (T_InMin > obj.m_paretoAggressiveness*T_OutMax)
                isParetoOptimal = 1;else, isParetoOptimal = 0; end
        end

        %% OUTPUT
        function obj = plotConvergence(obj)
            % Plot Pareto front
            obj.plotParetoFront();
        end

        function plotParetoFront(obj)
            plt = PlotId;
            figure(plt.pareto_front); set(gcf, 'Name', 'ParetoFront')
            plot(obj.m_history.constraint.volFrac(1:obj.m_iter),obj.m_history.objective(1:obj.m_iter),'-ko', 'LineWidth',2,'MarkerFaceColor','r');
            xlabel('Volume Fraction');
            if (strcmp(obj.m_objective,'compliance') == 1)
                ylabel('Normalized Compliance $C/C_0$');
            else
                ylabel('Normalized Objective $\varphi/\varphi_0$');
            end
            title('Pareto Front')
        end
    end


end

