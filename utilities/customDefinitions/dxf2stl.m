% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Code for exporting  .STL from a .DXF file                                 %
% %                                                                           %
% % This Matlab code was written by:                                          %
% % - Amir M. Mirzendehdel, Aerospace Engineering Department, KU              %
% % - Krishnan Suresh, Mechanical Engineering Department, UW-Madison          %
% %                                                                           %
% % Please send your comments to: amirzend@ku.edu                             %
% %                                                                           %
% % The code is intended for educational purposes and theoretical details     %
% % are discussed in the textbook:                                            %
% % Introduction to Shape and Topology Optimization using MATLAB              %
% %                                                                           %
% % Disclaimer:                                                               %
% % The authors reserves all rights but do not guaranty that the code is      %
% % free from errors. Furthermore, we shall not be liable in any event        %
% % caused by the use of the program.                                         %
% %                                                                           %
% % License:                                                                  %
% % This software is used, copied and distributed under the licensing         %
% % agreement contained in the file LICENSE in the top directory of           %
% % the distribution.                                                         %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function dxf2stl(dxf_filename,thickness,stl_filename)
% % dxf2stl
% % Convert .DXF file to .STL file
% % dxf_filename: input .dxf
% % thickness: input scalar
% % stl_filename: output .stl
% 
% if nargin < 2, thickness = 1; end
% 
% if nargin < 3
%     [path,example_name,~] = fileparts(dxf_filename);
%     stl_filename = [path example_name '.stl'];
% end
% %% 1. Import DXF
% dxf_filename
% [c_Line, ~, ~, ~, ~] = f_LectDxf(dxf_filename);
% rawLines = cell2mat(c_Line(:,1))
% 
% %% 2. Stitch Lines into Loops
% % Rounding to snap points together (adjust 4 to 2 if your scale is very small)
% V = round([rawLines(:,1:2); rawLines(:,4:5)], 4);
% [uniqueV, ~, idx] = unique(V, 'rows', 'stable');
% edges = reshape(idx, [], 2);
% 
% % Build a graph of all line connections
% G = graph(edges(:,1), edges(:,2));
% 
% % Find all connected components (each component is a potential loop)
% bins = conncomp(G);
% numLoops = max(bins);
% 
% % Initialize an empty polyshape
% pgon = polyshape();
% 
% for i = 1:numLoops
%     % Extract indices for this specific loop
%     nodeIdx = find(bins == i);
%     subG = subgraph(G, nodeIdx);
% 
%     % Order the nodes to form a path
%     loopPath = dfsearch(subG, 1);
% 
%     % Extract coordinates
%     x_loop = uniqueV(nodeIdx(loopPath), 1);
%     y_loop = uniqueV(nodeIdx(loopPath), 2);
% 
%     % Add this loop to our total shape
%     % polyshape automatically detects if a loop is inside another (a hole)
%     tempPgon = polyshape(x_loop, y_loop);
%     pgon = addboundary(pgon, tempPgon.Vertices);
% end
% 
% %% 3. Extrude to 3D
% tr2d = triangulation(pgon); % This triangulation respects the holes
% nodes2d = tr2d.Points;
% conn2d = tr2d.ConnectivityList;
% n = size(nodes2d, 1);
% 
% % Create 3D Vertices
% v3d = [nodes2d, zeros(n, 1); nodes2d, ones(n, 1) * thickness];
% 
% % Create Faces (Bottom and Top)
% fBottom = conn2d;
% fTop = conn2d + n;
% fTop = fTop(:, [1 3 2]); % Flip winding for outward normal
% 
% % Create Side Walls
% % freeBoundary on a triangulation with holes correctly identifies
% % both the outer edge and the inner hole edges.
% bEdges = freeBoundary(tr2d);
% fSides = zeros(size(bEdges,1)*2, 3);
% for i = 1:size(bEdges,1)
%     v1 = bEdges(i,1); v2 = bEdges(i,2);
%     v3 = v1 + n;      v4 = v2 + n;
%     fSides(2*i-1, :) = [v1, v2, v4];
%     fSides(2*i, :)   = [v1, v4, v3];
% end
% 
% %% 4. Export and Visualize
% TR = triangulation([fBottom; fTop; fSides], v3d);
% stlwrite(TR, stl_filename);
% 
% pltId = PlotId;
% fig = figure(pltId.geom_stl); clf(fig, 'reset');
% set(gcf, 'Name', 'STL');
% trisurf(TR, 'FaceColor', 'cyan');
% lighting gouraud
% material metal
% camlight
% axis equal; view(3); 
% end

function stl_filename = dxf2stl(dxf_filename, thickness, stl_filename)
% dxf2stl
% Convert .DXF file to a watertight extruded .STL
% Supports: LINE, LWPOLYLINE, POLYLINE/VERTEX
%
% dxf_filename: input .dxf
% thickness   : extrusion thickness (scalar)
% stl_filename: output .stl

if nargin < 2 || isempty(thickness), thickness = 1; end

if nargin < 3 || isempty(stl_filename)
    [p, name, ~] = fileparts(dxf_filename);
    stl_filename = fullfile(p, [name '.stl']);
end

%% 1) Read DXF into line segments (works for LINE + polyline)
rawLines = readDXFSegments_asLines6(dxf_filename);   % Nx6 [x1 y1 z1 x2 y2 z2]
if isempty(rawLines)
    error('No LINE/LWPOLYLINE/POLYLINE entities found (or could not parse).');
end

%% 2) Stitch lines into loops -> polyshape with holes
snapDigits = 4;  % increase if coordinates are noisy; decrease if very small scale
V = round([rawLines(:,1:2); rawLines(:,4:5)], snapDigits);
[uniqueV, ~, idx] = unique(V, 'rows', 'stable');
edges = reshape(idx, [], 2);

% Build graph and connected components
G = graph(edges(:,1), edges(:,2));
bins = conncomp(G);
numComps = max(bins);

pgon = polyshape();  % start empty

for ci = 1:numComps
    nodeIdx = find(bins == ci);

    % edges fully inside this component
    inComp = ismember(edges(:,1), nodeIdx) & ismember(edges(:,2), nodeIdx);
    subEdges = edges(inComp,:);

    % Trace an ordered loop (assumes this component is a single cycle)
    loopNodes = traceClosedLoop(nodeIdx, subEdges);

    x_loop = uniqueV(loopNodes, 1);
    y_loop = uniqueV(loopNodes, 2);

    % Add as boundary (polyshape will treat inside loops as holes automatically)
    if isempty(pgon.Vertices)
        pgon = polyshape(x_loop, y_loop);
    else
        pgon = addboundary(pgon, x_loop, y_loop);
    end
end

% Optional cleanup: merge tiny slivers etc.
pgon = simplify(pgon);

if isempty(pgon.Vertices)
    error('Failed to build a valid polyshape from DXF segments.');
end

%% 3) Triangulate 2D region and extrude
tr2d   = triangulation(pgon);    % respects holes
nodes2d = tr2d.Points;
conn2d  = tr2d.ConnectivityList;
n = size(nodes2d,1);

v3d = [nodes2d, zeros(n,1); nodes2d, thickness*ones(n,1)];

% bottom & top
fBottom = conn2d;
fTop    = conn2d + n;
fTop    = fTop(:, [1 3 2]); % flip winding for outward normal

% side walls from free boundary
bEdges = freeBoundary(tr2d); % boundary edges (outer + holes)
fSides = zeros(size(bEdges,1)*2, 3);
for i = 1:size(bEdges,1)
    v1 = bEdges(i,1); v2 = bEdges(i,2);
    v3 = v1 + n;      v4 = v2 + n;
    fSides(2*i-1, :) = [v1, v2, v4];
    fSides(2*i, :)   = [v1, v4, v3];
end

TR = triangulation([fBottom; fTop; fSides], v3d);

%% 4) Write STL + quick visualization
stlwrite(TR, stl_filename);

pltId = PlotId;
fig = figure(pltId.geom_stl); clf(fig, 'reset');
set(gcf, 'Name', 'STL Model');
trisurf(TR, 'FaceColor', 'cyan', 'EdgeColor','k');
axis equal; view(3);
camlight; lighting gouraud; material metal;

end

%% ======================= Helpers =======================

function rawLines6 = readDXFSegments_asLines6(filename)
% Returns Nx6 [x1 y1 z1 x2 y2 z2] from LINE/LWPOLYLINE/POLYLINE entities.

L = readlines(filename);
L = strip(L);
nL = numel(L);

segs = zeros(0,4);  % [x1 y1 x2 y2]

i = 1;
while i <= nL-1
    if L(i)=="0"
        ent = L(i+1);
        switch ent

            case "LINE"
                i = i + 2;
                x1 = NaN; y1 = NaN; x2 = NaN; y2 = NaN;
                while i <= nL-1 && L(i)~="0"
                    code = str2double(L(i));
                    val  = str2double(L(i+1));
                    switch code
                        case 10, x1 = val;
                        case 20, y1 = val;
                        case 11, x2 = val;
                        case 21, y2 = val;
                    end
                    i = i + 2;
                end
                if all(isfinite([x1 y1 x2 y2]))
                    segs(end+1,:) = [x1 y1 x2 y2]; %#ok<AGROW>
                end
                continue

            case "LWPOLYLINE"
                i = i + 2;
                V = zeros(0,2);
                pendingX = NaN;
                closed = false;

                while i <= nL-1 && L(i)~="0"
                    code = str2double(L(i));
                    valS = L(i+1);
                    i = i + 2;

                    switch code
                        case 70
                            flags  = str2double(valS);
                            closed = bitand(flags, 1) ~= 0;
                        case 10
                            pendingX = str2double(valS);
                        case 20
                            y = str2double(valS);
                            x = pendingX;
                            if isfinite(x) && isfinite(y)
                                V(end+1,:) = [x y]; %#ok<AGROW>
                            end
                            pendingX = NaN;
                        otherwise
                            % ignore bulge/width/etc
                    end
                end

                if size(V,1) >= 2
                    segs = [segs; [V(1:end-1,:) V(2:end,:)]]; %#ok<AGROW>
                    if closed
                        segs(end+1,:) = [V(end,:) V(1,:)]; %#ok<AGROW>
                    end
                end
                continue

            case "POLYLINE"
                % Old-style: POLYLINE then many VERTEX then SEQEND
                i = i + 2;

                closed = false;
                V = zeros(0,2);

                % read POLYLINE header until next 0
                while i <= nL-1 && L(i)~="0"
                    code = str2double(L(i));
                    valS = L(i+1);
                    if code == 70
                        flags  = str2double(valS);
                        closed = bitand(flags, 1) ~= 0;
                    end
                    i = i + 2;
                end

                % now expect 0 VERTEX blocks until 0 SEQEND
                while i <= nL-1
                    if L(i)=="0" && L(i+1)=="VERTEX"
                        i = i + 2;
                        x = NaN; y = NaN;
                        while i <= nL-1 && L(i)~="0"
                            code = str2double(L(i));
                            val  = str2double(L(i+1));
                            switch code
                                case 10, x = val;
                                case 20, y = val;
                            end
                            i = i + 2;
                        end
                        if all(isfinite([x y]))
                            V(end+1,:) = [x y]; %#ok<AGROW>
                        end
                    elseif L(i)=="0" && L(i+1)=="SEQEND"
                        i = i + 2;
                        break
                    else
                        i = i + 1;
                    end
                end

                if size(V,1) >= 2
                    segs = [segs; [V(1:end-1,:) V(2:end,:)]]; %#ok<AGROW>
                    if closed
                        segs(end+1,:) = [V(end,:) V(1,:)]; %#ok<AGROW>
                    end
                end
                continue

            otherwise
                % not an entity we care about
        end
    end

    i = i + 1;
end

% Convert segs Nx4 -> Nx6 with z=0
rawLines6 = [segs(:,1:2), zeros(size(segs,1),1), segs(:,3:4), zeros(size(segs,1),1)];
end

function loopNodes = traceClosedLoop(nodeIdx, subEdges)
% Try to order nodes around a single closed loop.
% Works best when every node has degree 2 (a cycle).

% Build local adjacency for all nodes in this component
maxNode = max(nodeIdx);
adj = cell(maxNode,1);
for k = 1:size(subEdges,1)
    a = subEdges(k,1); b = subEdges(k,2);
    adj{a}(end+1) = b; %#ok<AGROW>
    adj{b}(end+1) = a; %#ok<AGROW>
end

% pick a start node that exists
start = nodeIdx(1);

% If it's a clean cycle, degrees should be 2
deg = cellfun(@numel, adj(nodeIdx));
isCycle = all(deg == 2);

if isCycle
    loopNodes = zeros(numel(nodeIdx),1);
    loopNodes(1) = start;
    prev = 0;
    curr = start;

    for t = 2:numel(nodeIdx)
        nbrs = adj{curr};
        if prev == 0
            next = nbrs(1);
        else
            if nbrs(1) == prev
                next = nbrs(2);
            else
                next = nbrs(1);
            end
        end
        loopNodes(t) = next;
        prev = curr;
        curr = next;
    end
    % loopNodes currently ends at start's neighbor, not start; that's fine for polyshape
    return
end

% Fallback (non-cycle graph): return nodes as-is (may create a messy polygon)
warning('Component is not a simple cycle (node degree != 2 everywhere). Loop ordering may be incorrect.');
loopNodes = nodeIdx(:);
end
