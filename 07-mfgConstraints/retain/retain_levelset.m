%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discription:                                                              %
% This is a class for imposing edge retain constraints used in              %
% shape and topology optimization (e.g., level-set).                             %
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

classdef  retain_levelset < mfgConstraints
    properties (GetAccess = 'public', SetAccess = 'protected')
        m_retainElems; % elements to be retained
        m_retainNeighborIds; % for retain
    end
    % methods (Static,Access = private)
    %     function retainElems = dilateRetain(retainElems)
    %         se = strel('disk', 1);
    %         retainElems = imdilate(retainElems, se);
    %     end
    % end
    methods
        %% CONSTRUCTOR
        function obj = retain_levelset(solver,edges)
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
            % Makes the sensitivity of the retain element equal to zero
            % meaning the element would be kept
            filteredSensitivity = sensField;
            filteredSensitivity(obj.m_retainElems==1) = 0;

        end

        function obj = retainEdge(obj, seg)
            retainElems = zeros(obj.m_solver.m_ny, obj.m_solver.m_nx);

            nNode_y = obj.m_solver.m_ny + 1;

            for segId = 1:numel(seg)
                segNodes = obj.m_solver.findNodesOnEdge(seg(segId));

                for k = 1:numel(segNodes)
                    nodeId = segNodes(k);

                    [yNode, xNode] = obj.nodeId2sub(nodeId, nNode_y);

                    retainElems = obj.markElemsTouchingNode( ...
                        retainElems, yNode, xNode, obj.m_solver.m_existingElems);
                end
            end

            retainElems = obj.dilateRetain(retainElems);

            obj.m_retainElems = obj.m_retainElems | (retainElems == 1);
            obj.m_retainElems = obj.m_retainElems .* obj.m_solver.m_existingElems;
        end

        function obj = retainElements(obj)
            retainElems = zeros(obj.m_solver.m_ny, obj.m_solver.m_nx);
            nNode_y = obj.m_solver.m_ny + 1;

            for scenarioId = 1:obj.m_solver.m_numScenarios
                scenarioNodes = obj.m_solver.m_forcedNodes{scenarioId};

                for k = 1:numel(scenarioNodes)
                    nodeId = scenarioNodes(k);

                    % Correct nodeId -> (yNode,xNode)
                    [yNode, xNode] = obj.nodeId2sub(nodeId, nNode_y);

                    % Retain the up-to-4 elements touching this node
                    retainElems = obj.markElemsTouchingNode( ...
                        retainElems, yNode, xNode, obj.m_solver.m_existingElems);
                end
            end

            % Robust safety band (your helper uses disk radius 2 now)
            retainElems = obj.dilateRetain(retainElems);

            % Merge with existing retained set and mask by existingElems
            obj.m_retainElems = obj.m_retainElems | (retainElems == 1);
            obj.m_retainElems = obj.m_retainElems .* obj.m_solver.m_existingElems;
        end
    end
    methods (Static, Access = private)
        function [yNode, xNode] = nodeId2sub(nodeId, nNode_y)
            % nodeId = (xNode-1)*nNode_y + yNode  (MATLAB column-major)
            yNode = mod(nodeId-1, nNode_y) + 1;
            xNode = floor((nodeId-1)/nNode_y) + 1;
        end

        function retainElems = markElemsTouchingNode(retainElems, yNode, xNode, existingElems)
            % Elements touching node (yNode,xNode) are:
            % (yNode-1,xNode-1), (yNode-1,xNode), (yNode,xNode-1), (yNode,xNode)
            ny = size(existingElems,1);
            nx = size(existingElems,2);

            Js = [yNode-1, yNode-1, yNode,   yNode  ];
            Is = [xNode-1, xNode,   xNode-1, xNode  ];

            for t = 1:4
                j = Js(t); i = Is(t);
                if j>=1 && j<=ny && i>=1 && i<=nx && existingElems(j,i)==1
                    retainElems(j,i) = 1;
                end
            end
        end

        function retainElems = dilateRetain(retainElems)
            se = strel('disk', 2);
            retainElems = imdilate(retainElems, se);
        end
    end

end
