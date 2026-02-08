%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

classdef PlotBC
    properties
        %% elasticity
        fixed_U;
        fixed_V;
        force;
        acceleration;
        %% thermal
        fixed_T;
        flux;
        internal_heat;
        %% fluid
        noSlip_U = 'c';
        flow_U = [0.302, 0.745, 0.933];
        noSlip_V = 'm';
        flow_V = [0.9 0 0.9];
        fixed_P = [0.4940 0.1840 0.5560];

        activeDomain;
    end
    methods
        function obj = PlotBC()
            %% elasticity
            obj.fixed_U = struct();
            obj.fixed_U.color = 'k';
            obj.fixed_U.marker = 'x';

            obj.fixed_V = struct();
            obj.fixed_V.color = 'k';
            obj.fixed_V.marker = 'o';

            obj.force = struct();
            obj.force.color = [0 0.4470 0.7410];
            obj.force.marker = '^';

            obj.acceleration = struct();
            obj.acceleration.color = [0.4660 0.6740 0.1880];
            obj.acceleration.marker = '^';

            %% thermal
            obj.fixed_T = struct();
            obj.fixed_T.color = 'r';
            obj.fixed_T.marker = 'square';

            obj.flux = struct();
            obj.flux.color = [0.6350 0.0780 0.1840];
            obj.flux.marker = '^';

            obj.internal_heat = struct();
            obj.internal_heat.color = [0.9290 0.6940 0.1250];
            obj.internal_heat.marker = '^';

            %% fluid
            obj.noSlip_U = struct();
            obj.noSlip_U.color = 'c';
            obj.noSlip_U.marker = 'o';

            obj.flow_U = struct();
            obj.flow_U.color = [0.302, 0.745, 0.933];
            obj.flow_U.marker = '>';

            obj.noSlip_V = struct();
            obj.noSlip_V.color = 'm';
            obj.noSlip_V.marker = '*';

            obj.flow_V = struct();
            obj.flow_V.color = [0.9 0 0.9];
            obj.flow_V.marker = '^';

            obj.fixed_P = struct();
            obj.fixed_P.color = [0.4940 0.1840 0.5560];
            obj.fixed_P.marker = 'diamond';  

            obj.activeDomain = struct();
            obj.activeDomain.color = [0.8 0.8 0.8];
            obj.activeDomain.marker = 'square';
        end
    end
end



