%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discription:                                                              %
% This is a class for imposing edge retain constraints used in              %
% shape and topology optimization via topological sensitivity field (TSF).  %
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

classdef  retain_tsf < mfgConstraints
    properties (GetAccess = 'public', SetAccess = 'protected')
        m_retainElems; % elements to be retained
        m_retainNeighborIds; % for retain
    end
    methods (Static,Access = private)
        function retainElems = dilateRetain(retainElems)
            se = strel('disk', 1);
            retainElems = imdilate(retainElems, se);
        end
    end
    methods
        %% CONSTRUCTOR
        function obj = retain_tsf(solver,edges)
            % check if solver is valid
            if (~isa(solver, 'simulation2d')), error('solver must be an instance of simulation2d class!');end % check if solver is valid

            % constructor based on superclass
            obj = obj@mfgConstraints(solver);

            if (isempty(edges))
                error('edges array is empty!');
            end

            % initialize
            obj.m_retainElems = zeros(size(obj.m_solver.m_existingElems));

            % retain neighborhood
            obj.m_retainNeighborIds = [
                0  0;  % Center (current element)
                0 -1;  % Left
                -1  0; % Below
                0  1;  % Right
                1  0;  % Above
                ];

            obj = obj.retainEdge(edges);
        end

        function [filteredDesign] = filterDesign(obj, design)
            % filter the design variables
            % input: obj, design variables
            % output: obj, filtered design variables

            % Apply retain on density
            % Makes the density of the retain element equal to 1
            filteredDesign = design;
            filteredDesign(obj.m_retainElems==1) = 1;

        end

        function [filteredSensitivity] = filterSensitivity(obj, ~, sensField)
            % filter the sensitivity fields
            % input: obj, sensitivity fields
            % output: obj, filtered sensitivity fields

            % Apply retain on sensitivity
            % Makes the sensitivity of the retain element equal to the minimum across all elements
            % meaning the element would be kept
            filteredSensitivity = sensField;
            filteredSensitivity(obj.m_retainElems==1) = max(sensField(:));

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

        function obj = retainElements(obj)
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
    end
end
