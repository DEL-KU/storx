%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Discription:                                                              %
% This code defines a class named 'brep2d' that represents a 2D boundary    %
% representation (b-rep) model. It includes properties for storing various  %
% attributes of the b-rep, such as vertices, segments, area, perimeter,     %
% and bounding box dimensions. The class provides methods for reading b-rep %
% data from a file, calculating the area under line segments, and           %
% converting the b-rep to a format suitable for solving partial             %
% differential equations (PDEs) in MATLAB.                                  %
% The constructor initializes the object by either reading from a file or   %
% using a provided b-rep structure, computes the polygon representation,    %
% and calculates the model's center, scale, bounding box, area, and         %
% perimeter.                                                                %
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

classdef brep2d < handle
    properties(GetAccess = 'public', SetAccess = 'public')
        % m_brep is a struct with the following required fields:
        %   - vertices  : 2 x N double array of vertex coordinates
        %   - segments  : 4 x M double array of segment data
        m_brep;

        % m_pdeGeom is a struct or other data representation used for
        %    creating the PDE Toolbox geometry and plotting (e.g. `decsg`)
        m_pdeGeom;

        % m_pdeModel is a PDE model object created using the PDE Toolbox
        %    (e.g. `createpde`), which is used for solving PDEs on the geometry.
        %    It may contain fields like:
        %       - geometry : the geometry representation used in the PDE model
        %       - mesh     : the mesh generated for the geometry
        %       - boundary  : boundary conditions applied to the geometry
        %       - coefficients : coefficients for the PDE equations
        %       - solution  : solution to the PDE problem defined on the geometry
        %    This field is typically used for solving PDEs on the B-rep geometry.
        %    It may be empty if the B-rep is not yet associated with a PDE model.
        %    Example: m_pdeModel = createpde(obj.m_pdeGeom);
        m_pdeModel;

        % m_polygon is a structure representing the polygon corresponding
        %    to this B-rep, typically used for point-in-polygon queries.
        %    It may contain fields like:
        %       - vertices : 2 x N double array of polygon vertices
        %       - edges    : indices of polygon edges
        m_polygon;

        % m_numSegsPerArc is an integer specifying the number of straight
        %    subdivisions used to approximate each arc segment.
        m_numSegsPerArc;

        % m_brepSegsMapping is an index array mapping each B-rep segment
        %    to its corresponding PDE geometry segment.
        m_brepSegsMapping;

        % m_numBndrySegs is the number of boundary segments comprising the B-rep.
        m_numBndrySegs;

        % m_segLengths is a 1 x M double array containing the length of each boundary segment.
        m_segLengths;

        % m_modelScale is a scalar representing the characteristic scale
        %    of the model, e.g. average bounding-box size.
        m_modelScale;

        % m_boundingBox is a 2 x 2 double array representing the min/max
        %    corner of the B-rep bounding box:
        %       [ xmin  xmax ;
        %         ymin  ymax ]
        m_boundingBox;

        % m_boxSizes is a 1 x 3 double array containing:
        %    [ boxWidth  boxHeight  1 ] to describe box size in x and y
        m_boxSizes;

        % m_area is a scalar containing the computed area enclosed by the B-rep.
        m_area;

        % m_perimeter is a scalar containing the total perimeter length of the B-rep.
        m_perimeter;

        % m_center is a 2 x 1 vector containing the geometric center of the B-rep
        %    computed as the mean of its vertex coordinates.
        m_center;
    end


    methods(Static)
        %% READ B-REP FILE
        function brep = readBrep2D(fileName)
            % Read a 2D boundary representation (b-rep) from a file
            fid = fopen(fileName,'r'); % Open the file for reading
            nPoints = fscanf(fid,' %d ',1); % Read the number of points
            brep.vertices = fscanf(fid,' %f ',[2, nPoints]); % Read the vertex coordinates
            nSegments = fscanf(fid,' %d ',1); % Read the number of segments
            brep.segments = fscanf(fid,' %f ',[4, nSegments]); % Read the segment data
            fclose(fid); % Close the file
        end
        function area = areaUnderLineSegment(startPt,endPt,y_min)
            if(nargin < 3), y_min = 0; end % Default y_min to 0 if not provided
            % Calculate the area under a line segment defined by two points
            base = (startPt(1)-endPt(1)); % Calculate the base length
            a1 = base*(startPt(2)-y_min); % Area of the rectangle formed by the start point
            a2 = 0.5*(base)*(endPt(2)-startPt(2)); % Area of the triangle formed by the end point
            area = a1 + a2; % Total area under the line segment
        end
        function validateBrep(brep)
            % VALIDATEBREP Function to validate the properties of a boundary representation (BRep)
            %
            % Input Arguments:
            %     brep - structure containing the BRep properties to validate
            %
            % Throws:
            %     Assertion errors if the properties do not meet the specified criteria

            % Check if 'vertices' field exists and is a 2 x N numeric array
            assert(isfield(brep, 'vertices') && isnumeric(brep.vertices) && size(brep.vertices,1)==2, ...
                'm_brep.vertices must be a 2 x N numeric array.');
            % Check if 'segments' field exists and is a 4 x M numeric array
            assert(isfield(brep, 'segments') && isnumeric(brep.segments) && size(brep.segments,1)==4, ...
                'm_brep.segments must be a 4 x M numeric array.');
        end
    end

    methods
        %% CONSTRUCTOR
        function obj = brep2d(brep)
            % BREP2D Constructor for the brep2d class
            % This constructor initializes the brep2d object with a boundary representation (B-rep).
            % Input Arguments:
            %     brep - a structure containing the B-rep data or a file name
            %            % The B-rep structure should have the following fields:
            %   - vertices: a 2 x N double array of vertex coordinates
            %   - segments: a 4 x M double array of segment data
            %            % The segments array should have the following format:
            %   - segments(1,:) : segment type (1 for line, 2 for arc, -1 for virtual segments)
            %   - segments(2,:) : start vertex index
            %   - segments(3,:) : end vertex index
            %   - segments(4,:) : additional data (e.g., center vertex for arcs)
            %            % Example:
            %   brep.vertices = [x1, x2, ..., xN; y1, y2, ..., yN];
            %   brep.segments = [type1, type2, ..., typeM;
            %                            start1, start2, ..., startM;
            %                            end1, end2, ..., endM;
            %                            center1, center2, ..., centerM]; % for arcs only: positive for clockwise and negative for counter-clockwise

            if (ischar(brep)), obj.m_brep = obj.readBrep2D(brep);
            else, obj.m_brep = brep; end

            % Validate the boundary representation (B-rep) to ensure it meets the required criteria
            obj.validateBrep(obj.m_brep);

            % Store the number of boundary segments and compute the polygon representation
            obj.m_numBndrySegs = size(obj.m_brep.segments,2);

            obj.m_numSegsPerArc = 50; % Set the number of segments per arc for approximation

            % Calculate the center and scale of the model based on vertex coordinates
            x = obj.m_brep.vertices(1,:); y = obj.m_brep.vertices(2,:);
            obj.m_center = [mean(x);mean(y)];
            obj.m_modelScale = ((max(x)-min(x))+(max(y)-min(y)))/2;

            % Define the bounding box of the model
            obj.m_boundingBox = [min(x) max(x); min(y) max(y)];
            obj.m_segLengths = zeros(1,obj.m_numBndrySegs); % Initialize segment lengths

            % Convert BREP to PDE geometry representation
            obj = obj.convertBreptoPdeGeom();

            % Compute the polygon representation of the B-rep
            obj = obj.convertPdeGeomtoPolygon();

            % Calculate the sizes of the bounding box
            lx = obj.m_boundingBox(1,2) - obj.m_boundingBox(1,1);
            ly = obj.m_boundingBox(2,2) - obj.m_boundingBox(2,1);
            obj.m_boxSizes = [lx ly 1 ]; % Set the box sizes

            % Compute the area and perimeter of the BREP
            obj.m_area = obj.brepArea();
            obj.m_perimeter = sum(obj.m_segLengths);
        end

        %% CONVERT B-REP TO MATLAB PDEGEOM
        function obj =  convertBreptoPdeGeom(obj)
            % CONVERTBREPTOPDEGEOM Convert from brep format to Matlab's pdegeom format
            %
            % Input Arguments:
            %     obj - object containing the brep data
            %
            % Output Arguments:
            %     obj - updated object with pdegeom format

            brep = obj.m_brep; % Extract the brep data from the object
            v = brep.vertices; % Get the vertices from the brep
            nSegments = size(brep.segments,2); % Number of segments in the brep
            pdeGeom = zeros(11,1); % Initialize pdeGeom array
            segmentMapping = []; % Initialize segment mapping array
            npdeGeomSegments = 0; % Counter for pdeGeom segments
            dl = 1; dr = 0; % Define left and right segment lengths

            % Loop through each segment in the brep
            for i = 1:nSegments
                breptype = brep.segments(1,i); % Get the type of the segment
                vs = brep.segments(2,i); ve = brep.segments(3,i); % Start and end vertices
                xs = v(1,vs); xe = v(1,ve); % X-coordinates of start and end
                ys = v(2,vs); ye = v(2,ve); % Y-coordinates of start and end

                if (breptype == 1) % line
                    pdetype = 2; % Set pde type for line
                    npdeGeomSegments = npdeGeomSegments+1; % Increment segment count
                    pdeGeom(:,npdeGeomSegments) = ...
                        [pdetype xs xe ys ye dl dr 0 0 0 0]; % Store line segment data
                    segmentMapping(npdeGeomSegments) = i; %#ok
                    obj.m_segLengths(i) = sqrt( (xe-xs)^2 + (ye-ys)^2); % Calculate segment length

                elseif (breptype == 2) % arc
                    segtype = 1;                 % 1 = circle segment (dl format)
                    vc = abs(brep.segments(4,i));
                    xc = v(1,vc);  yc = v(2,vc);
                    R  = hypot(xc - xs, yc - ys);
                    d = sqrt((xe-xs)^2+(ye-ys)^2); % Arc width
                    % Ensure CCW: if your B-rep says it's clockwise, reverse direction
                    if brep.segments(4,i) > 0    % (per your convention) clockwise
                        % reverse endpoints, and swap left/right because direction flipped
                        x1 = xe; y1 = ye;
                        x2 = xs; y2 = ys;
                        left  = dr;
                        right = dl;
                    else                         % already CCW
                        x1 = xs; y1 = ys;
                        x2 = xe; y2 = ye;
                        left  = dl;
                        right = dr;
                    end

                    npdeGeomSegments = npdeGeomSegments + 1;
                    segmentMapping(npdeGeomSegments) = i; %#ok

                    pdeGeom(:,npdeGeomSegments) = [ ...
                        segtype;  x1; x2;  y1; y2;  left; right;  xc; yc;  R; R ];

                    % Calculate the angle of the arc
                    theta = real(2*asin(d/2/R));
                    arcLength = R*theta; % Calculate arc length
                    obj.m_segLengths(i) = arcLength; % Store arc length
                    numArcSegments = obj.m_numSegsPerArc; % Number of segments per arc
                    vecStart = [xs-xc ys-yc]; % Vector from center to start
                    vecStart = vecStart/norm(vecStart); % Normalize vector
                    vecEnd = [xe-xc ye-yc]; % Vector from center to end
                    vecEnd = vecEnd/norm(vecEnd); % Normalize vector
                    thetaStart = atan2(vecStart(2),vecStart(1)); % Start angle
                    thetaEnd = atan2(vecEnd(2),vecEnd(1)); % End angle

                    % Adjust angles based on direction of the arc
                    if (brep.segments(4,i) > 0) % clock-wise
                        if (thetaEnd > thetaStart) % not allowed
                            thetaEnd = thetaEnd - 2*pi; % Adjust end angle
                        end
                    else % counter clock-side
                        if (thetaEnd < thetaStart) % not allowed
                            thetaEnd = thetaEnd + 2*pi; % Adjust end angle
                        end
                    end

                    % Loop through the number of arc segments
                    for j = 1:numArcSegments
                        % npdeGeomSegments = npdeGeomSegments+1; % Increment segment count
                        % segmentMapping(npdeGeomSegments) = i; %#ok
                        theta0 = thetaStart + ...
                            (j-1)*(thetaEnd-thetaStart)/numArcSegments; % Start angle for segment
                        
                        xs = xc + R*cos(theta0); % Calculate segment endpoints
                        ys = yc + R*sin(theta0);

                        % Update bounding box
                        obj.m_boundingBox(1,1) = min(obj.m_boundingBox(1,1),xs);
                        obj.m_boundingBox(1,2) = max(obj.m_boundingBox(1,2),xs);
                        obj.m_boundingBox(2,1) = min(obj.m_boundingBox(2,1),ys);
                        obj.m_boundingBox(2,2) = max(obj.m_boundingBox(2,2),ys);
                    end

                elseif (breptype == 3) % circle: same start and end vertices
                    pdetype = 1; % Set pde type for arc
                    vc = abs(brep.segments(4,i)); % Get vertex center for circle
                    xc = v(1,vc); yc  = v(2,vc); % Circle center coordinates
                    R = sqrt((xc-xs)^2+(yc-ys)^2); % Circle radius
                    numArcSegments = obj.m_numSegsPerArc; % Number of segments for circle

                    % Loop through the number of arc segments for the circle
                    for j = 1:numArcSegments
                        npdeGeomSegments = npdeGeomSegments+1; % Increment segment count
                        segmentMapping(npdeGeomSegments) = i; %#ok
                        thetaStart = (j-1)*2*pi/numArcSegments; % Start angle for segment
                        thetaEnd = j*2*pi/numArcSegments; % End angle for segment
                        xs = xc + R*cos(thetaStart); xe = xc + R*cos(thetaEnd); % Calculate segment endpoints
                        ys = yc + R*sin(thetaStart); ye = yc + R*sin(thetaEnd);

                        % Update bounding box
                        obj.m_boundingBox(1,1) = min(obj.m_boundingBox(1,1),xs);
                        obj.m_boundingBox(1,2) = max(obj.m_boundingBox(1,2),xs);
                        obj.m_boundingBox(2,1) = min(obj.m_boundingBox(2,1),ys);
                        obj.m_boundingBox(2,2) = max(obj.m_boundingBox(2,2),ys);

                        % Store circle segment data
                        pdeGeom(:,npdeGeomSegments) = ...
                            [pdetype xs xe ys ye dl dr xc yc R R];
                    end
                end
            end

            obj.m_brepSegsMapping = segmentMapping; % Store the segment mapping
            
            % Clean-up
            [obj.m_pdeGeom, ~] = pdegeom_fixArcsProjectToCircle(pdeGeom);
            obj.m_pdeGeom = pdegeom_normalizeTo10Rows(obj.m_pdeGeom);

            % Convert pdeGeom to pdeModel
            % Create a PDE model using the geometry defined in pdeGeom
            obj.m_pdeModel = createpde(1);
            geometryFromEdges(obj.m_pdeModel,obj.m_pdeGeom);
        end

        function obj = convertPdeGeomtoPolygon(obj)
            % CONVERTPDEGEOMTOPOLYGON Convert from pdeGeom format to polygon representation
            % Input Arguments:
            %     obj - object containing the pdeGeom data
            % Output Arguments:
            %     obj - updated object with polygon representation
            % Extract vertices and segments from the pdeGeom

            brep = obj.m_brep;
            v = brep.vertices;
            nSegments = size(brep.segments,2);
            xv(1,1) = v(1,1);
            yv(1,1) = v(2,1);
            counter = 2;
            for i = 1:nSegments
                if(brep.segments(1,i)== -1)
                    xv(counter) = NaN;
                    yv(counter) = NaN;
                    counter = counter + 1;
                    node = brep.segments(3,i);
                    xv(counter) = v(1,node);
                    yv(counter) = v(2,node);
                    counter = counter+1;
                elseif (brep.segments(1,i)== 1)
                    node = brep.segments(3,i);
                    xv(counter) = v(1,node);
                    yv(counter) = v(2,node);
                    counter = counter+1;
                elseif (brep.segments(1,i)== 2)

                    node1 = brep.segments(2,i);
                    node2 = brep.segments(3,i);

                    center = abs(brep.segments(4,i));
                    dir = -sign(brep.segments(4,i));
                    [X,Y] = obj.PtsOnArc(node1,node2,center,obj.m_numSegsPerArc,dir);
                    xv(counter:counter+obj.m_numSegsPerArc-1) = X;
                    yv(counter:counter+obj.m_numSegsPerArc-1) = Y;
                    counter = counter + obj.m_numSegsPerArc;
                end

            end
            polygon0 = [xv;yv];
            polygon = zeros(size(polygon0));
            k = 0;
            counter = 0;
            % outside boundary
            for i = 1:size(polygon,2)
                if (isnan(polygon0(:,i)))
                    k = k+1;
                    continue;
                end
                if (mod(k,2)==0)
                    counter = counter + 1;
                    polygon(:,counter) = polygon0(:,i);
                end
            end
            % holes
            k = 0;
            for i = 1:size(polygon,2)
                if (isnan(polygon0(:,i)))
                    k = k+1;
                    counter = counter + 1;
                    polygon(:,counter) = polygon0(:,i);
                    continue;
                end
                if (mod(k,2)==1)
                    counter = counter + 1;
                    polygon(:,counter) = polygon0(:,i);
                end
            end
            % fix
            if (isnan(polygon(:,size(polygon,2))))
                polygon(:,size(polygon,2)) = [];
            end

            obj.m_polygon = [polygon polygon(:,1)]; % Close polygon()
        end


        %% SAMPLE POINTS ON ARC
        function [X,Y] = PtsOnArc(obj,node1,node2,center,nSections,dir)
            % Calculate points on an arc between two nodes around a center
            v = obj.m_brep.vertices;

            % Extract coordinates
            x1 = v(1,node1);    y1 = v(2,node1);
            x2 = v(1,node2);    y2 = v(2,node2);
            xc = v(1,center);   yc = v(2,center);

            % Radius vectors
            r1 = [x1 - xc; y1 - yc];
            r2 = [x2 - xc; y2 - yc];
            radius = norm(r1);

            % Compute initial angles using atan2 for robustness
            startAngle = atan2(r1(2), r1(1));
            endAngle = atan2(r2(2), r2(1));

            % Ensure angles are in the correct range and direction
            angleDiff = endAngle - startAngle;

            if dir == 1  % Counter-clockwise
                if angleDiff <= 0
                    angleDiff = angleDiff + 2*pi;
                end
            elseif dir == -1  % Clockwise
                if angleDiff >= 0
                    angleDiff = angleDiff - 2*pi;
                end
            else
                error('dir must be either 1 (CCW) or -1 (CW)');
            end

            % Generate angles along arc
            angles = linspace(startAngle, startAngle + angleDiff, nSections + 1);

            % Calculate points on arc
            X = xc + radius * cos(angles);
            Y = yc + radius * sin(angles);

            % Remove the first point (starting node already known)
            X = X(2:end);
            Y = Y(2:end);
        end

        function [sd,closestPts,closestSegs] = signedDistance(obj,pts)
            % SIGNEDDISTANCE Computes the signed distance from points to a boundary representation (BRep)
            %
            % Input Arguments:
            %     obj  - Object containing BRep and distance methods
            %     pts  - 2D points for which the signed distance is calculated
            %
            % Output Arguments:
            %     sd          - Signed distances of points to the BRep
            %     closestPts  - Closest points on the BRep to the input points
            %     closestSegs  - Segments of the BRep corresponding to the closest points

            % Determine if points are inside the BRep
            in = obj.inBrep(pts);
            % Initialize signed distance array
            sd = ones(size(in));
            % Assign negative distance for points outside the BRep
            sd(in==0) = -1;
            % Calculate minimum distances and closest points/segments
            [dMin,closestPts,closestSegs] = obj.distOfPointsToBrep(pts);
            % Update signed distances with minimum distances
            sd = sd .* dMin;
        end

        %% FIND POINTS IN THE POLYGON
        function in = inBrep(obj,pts)
            % INBREP Determines point membership within the polygon using inpolygon
            %
            % Input Arguments:
            %     obj  - Object containing polygon data
            %     pts  - 2D points to check for membership
            %
            % Output Arguments:
            %     in   - Logical array indicating if points are inside (1) or outside (0) the polygon

            % uses inpolygon to compute point membership
            in = inpolygon(pts(1,:),pts(2,:),obj.m_polygon(1,:),obj.m_polygon(2,:));
        end
        %% FIND DISTANCE OF POINT TO B-REP (NEEDED FOR MESHING)
        function [dMin,closestPts,closestSegs] = distOfPointsToBrep(obj,pts,segmentsToCheck)
            % DISTOFPOINTSTOBREP Function to calculate the minimum distance from points to a boundary representation (BRep)
            %
            % Input Arguments:
            %     obj - object containing BRep data
            %     pts - matrix of points to check distances to the BRep
            %     segmentsToCheck - optional array of segment indices to check
            %
            % Output Arguments:
            %     dMin - minimum distances from points to the BRep
            %     closestPts - closest points on the BRep to the input points
            %     closestSegs - indices of the closest segments in the BRep

            brep = obj.m_brep; % Access the BRep data from the object
            nSegments = size(brep.segments,2); % Get the number of segments in the BRep
            if (nargin < 3)
                segmentsToCheck = 1:nSegments; % Default to check all segments if none specified
            end
            nPts = size(pts,2); % Get the number of points
            dMin = 1e12*ones(1,nPts); % Initialize minimum distances to a large value
            closestPts = pts; % Default closest points are the input points
            closestSegs = zeros(1,nPts); % Initialize closest segment indices to zero
            for seg = segmentsToCheck
                breptype = brep.segments(1,seg); % Determine the type of segment
                if (breptype == 1) % line
                    [dSeg,closestPtsSeg] = obj.distOfPointsToLineSegment(pts,seg); % Calculate distances for line segments
                elseif (breptype == 2) % arc
                    [dSeg,closestPtsSeg] = obj.distOfPointsToArcSegment(pts,seg); % Calculate distances for arc segments
                end
                for n = 1:nPts
                    if (dSeg(n) < dMin(n)) % Check if the current segment provides a closer point
                        dMin(n) = dSeg(n); % Update minimum distance
                        closestPts(:,n) = closestPtsSeg(:,n); % Update closest point
                        closestSegs(n) = seg; % Update closest segment index
                    end
                end
            end
        end
        %% FIND DISTANCE OF POINT TO LINE SEGMENT
        function [d,closestPts] = distOfPointsToLineSegment(obj,pts,seg)
            % DISTOFPOINTSTOLINESEGMENT Calculate the distance from points to a line segment
            %
            % Input Arguments:
            %     obj   - object containing the BREP data
            %     pts   - matrix of points to measure distance from
            %     seg   - segment index in the BREP
            %
            % Output Arguments:
            %     d          - distances from points to the line segment
            %     closestPts - closest points on the line segment to the input points

            % Extract vertices and line segment endpoints
            brep = obj.m_brep;  v = brep.vertices;
            lineStart = v(:,brep.segments(2,seg));
            lineEnd = v(:,brep.segments(3,seg));
            ABS_TOL = 1e-10; % Define a small tolerance for length comparison
            L = norm(lineStart-lineEnd); % Calculate the length of the line segment
            closestPts = pts; % Initialize closest points to input points
            if (L < ABS_TOL) % Check if the line segment is effectively a point
                d = norm(pts-lineStart); % Distance to the point
                closestPts(1,:) = lineStart(1); % Closest point x-coordinate
                closestPts(2,:) = lineStart(2); % Closest point y-coordinate
                return; % Exit if the segment is a point
            end
            % Calculate the unit tangent vector of the line segment
            lineTangent = (lineEnd-lineStart)/L;
            v1(1,:) = pts(1,:) - lineStart(1); % Vector from line start to points
            v1(2,:) = pts(2,:) - lineStart(2);
            v2(1,:) = pts(1,:) - lineEnd(1); % Vector from line end to points
            v2(2,:) = pts(2,:) - lineEnd(2);
            numer = (lineEnd-lineStart)'*v1; % Numerator for projection
            u = numer/L^2; % Projection factor
            distToStart = norm(v1); % Distance from points to line start
            distToEnd= norm(v2); % Distance from points to line end
            distToLineSeg = abs((-v1(1,:)*lineTangent(2) + v1(2,:)*lineTangent(1))); % Perpendicular distance to line segment
            % Calculate distances based on projection factor
            d = (u < 0).*distToStart + (u >1).*distToEnd + (u>=0).*(u<= 1).*distToLineSeg;
            % Determine closest points on the line segment
            closestPts(1,:) = (u < 0).*lineStart(1) + (u >1).*lineEnd(1) + ...
                (u>=0).*(u<= 1).*((1-u).*lineStart(1) + u.*lineEnd(1));
            closestPts(2,:) = (u < 0).*lineStart(2) + (u >1).*lineEnd(2) + ...
                (u>=0).*(u<= 1).*((1-u).*lineStart(2) + u.*lineEnd(2));
        end
        %% FIND DISTANCE OF POINT TO ARC SEGMENT
        function [d,closestPts] = distOfPointsToArcSegment(obj,pts,seg)
            % DISTOFPOINTSTOARCSEGMENT Calculate the distance from points to an arc segment
            %
            % Input Arguments:
            %     obj   - object containing the BREP data
            %     pts   - matrix of points to measure distance from
            %     seg   - segment index in the BREP
            %
            % Output Arguments:
            %     d          - distances from points to the arc segment
            %     closestPts - closest points on the arc segment to the input points

            % Extract BREP data and arc segment parameters
            brep = obj.m_brep;
            vc = abs(brep.segments(4,seg)); % Vertex index for the center of the arc
            v = brep.vertices; % Vertex coordinates
            xc = v(1,vc); yc  = v(2,vc); % Center coordinates of the arc
            vs = brep.segments(2,seg); % Start vertex index
            ve = brep.segments(3,seg); % End vertex index
            xs = v(1,vs); ys = v(2,vs); % Start coordinates
            xe = v(1,ve); ye = v(2,ve); % End coordinates
            R = sqrt((xc-xs)^2+(yc-ys)^2); % Radius of the arc
            % Normalize vectors from center to start and end points
            vecStart = [xs-xc ys-yc];
            vecStart = vecStart/norm(vecStart);
            vecEnd = [xe-xc ye-yc];
            vecEnd = vecEnd/norm(vecEnd);
            % Calculate angles for the start and end of the arc
            thetaStart = atan2(vecStart(2),vecStart(1));
            thetaEnd = atan2(vecEnd(2),vecEnd(1));
            % Adjust angles based on the direction of the arc
            if (brep.segments(4,seg) > 0) % Clockwise
                if (thetaEnd > thetaStart) % Adjust if not allowed
                    thetaEnd = thetaEnd - 2*pi;
                end
            else % Counterclockwise
                if (thetaEnd < thetaStart) % Adjust if not allowed
                    thetaEnd = thetaEnd + 2*pi;
                end
            end
            % Sample the arc
            N  = 2*obj.m_numSegsPerArc; % Number of sample points
            theta = thetaStart + (0:N)*(thetaEnd-thetaStart)/N; % Angle samples
            arcPts = zeros(2,N+1); % Preallocate arc points
            arcPts(1,:) = xc + R*cos(theta); % X-coordinates of arc points
            arcPts(2,:)  = yc + R*sin(theta); % Y-coordinates of arc points
            % Initialize closest points and distances
            closestPts = pts;
            d = zeros(1,size(pts,2));
            % Find closest point on the arc for each input point
            for i = 1:size(pts,2)
                dArc = sqrt( (pts(1,i) - arcPts(1,:)).^2 +  (pts(2,i) - arcPts(2,:)).^2); % Distances to arc points
                [d(i),index] = min(dArc); % Minimum distance and index of closest point
                closestPts(:,i) = arcPts(:,index); % Store closest point
            end
        end

        function area = areaUnderArcSegment(obj, startPt, endPt, centerPt,dir, y_min,nSamples)
            % Approximate area under an arc by fine line segmentation.
            % nSamples = number of subdivisions along the arc.
            if nargin < 6, y_min = 0; end
            if nargin < 7, nSamples = 20; end  % resolution

            % Vectors from center to endpoints
            v1 = startPt(:) - centerPt(:);
            v2 = endPt(:)   - centerPt(:);
            R  = 0.5*(norm(v1)+norm(v2));

            % Start/end angles
            t1 = atan2(v1(2), v1(1));
            t2 = atan2(v2(2), v2(1));

            % Adjust angles based on the direction of the arc
            if (dir > 0) % Clockwise
                if (t2 > t1) % Adjust if not allowed
                    t2 = t2 - 2*pi;
                end
            else % Counterclockwise
                if (t2 < t1) % Adjust if not allowed
                    t2 = t2 + 2*pi;
                end
            end

            % Parametric subdivision
            ts = linspace(t1, t2, nSamples+1);
            pts = centerPt(:).' + R*[cos(ts(:)), sin(ts(:))];

            % Sum line contributions
            area = 0;
            for k = 1:nSamples
                area = area + obj.areaUnderLineSegment(pts(k,:), pts(k+1,:), y_min);
            end
        end

        function area = brepArea(obj)
            % Calculate the area of the BREP
            % This function computes the area enclosed by the BREP segments
            % Input Arguments:
            %     obj - object containing the BREP data
            %
            % Output Arguments:
            %     area - area enclosed by the BREP segments

            brep = obj.m_brep;
            v = brep.vertices;
            y_min = min(obj.m_boundingBox(2,:));
            nSegments = size(brep.segments,2);
            area = 0;
            for i = 1:nSegments
                breptype = brep.segments(1,i);
                vs = brep.segments(2,i);
                ve = brep.segments(3,i);
                if (breptype == 1) % line
                    area = area + obj.areaUnderLineSegment(v(:,vs),v(:,ve),y_min);
                elseif (breptype == 2) % arc
                    vc = abs(brep.segments(4,i));
                    dir = sign(brep.segments(4,i));
                    a = obj.areaUnderArcSegment(v(:,vs),v(:,ve),v(:,vc),dir,y_min);
                    area = area + a;

                elseif (breptype == 3) % circle
                    vc = abs(brep.segments(4,i));
                    radius = norm(v(:,vs)-v(:,vc));
                    if (i == 1)
                        area = area + pi*radius^2; %external boundary
                    else
                        area = area - pi*radius^2; %internal hole
                    end
                elseif (breptype == -1) % internal segment
                    % do nothing
                else
                    disp('Unknown breptype type');
                end
            end
        end
        %% NORMAL
        function normalVec = normalOfSegment(obj,seg)
            % NORMALOFSEGMENT Function to calculate the normal vector of a line segment
            %
            % Input Arguments:
            %     obj - object containing the boundary representation (brep)
            %     seg - index of the segment for which the normal vector is calculated
            %
            % Output Arguments:
            %     normalVec - calculated normal vector of the segment

            % Access the boundary representation and vertices
            brep = obj.m_brep;  v = brep.vertices;
            % Get the start and end points of the line segment
            lineStart = v(:,brep.segments(2,seg));
            lineEnd = v(:,brep.segments(3,seg));
            % Calculate the length of the line segment
            L = norm(lineStart-lineEnd);
            % Compute the unit tangent vector of the line segment
            lineTangent = (lineEnd-lineStart)/L;
            % Calculate the normal vector by rotating the tangent vector
            normalVec = [-lineTangent(2);lineTangent(1)];
        end
        %% PLOT B-REP WITH/WITHOUT LABELS
        function plotGeometryWithLabels(obj,figId)
            % PLOTGEOMETRYWITHLABELS Function to plot geometry with labels
            %
            % Input Arguments:
            %     obj   - object containing geometry data
            %     figId - figure ID for plotting

            plt = PlotId;
            switch (nargin)
                case 1
                    % Default to brep figure ID if only object is provided
                    figId = plt.brep;
                otherwise
                    % use input values
            end
            % Call the plotGeometry function to perform the actual plotting
            obj.plotGeometry(figId,1);
        end
        function plotGeometry(obj,figId,plotLabels,title)
            % Plot geometry with exactly one central label per unique brep segment

            plt = PlotId;
            if nargin < 2, figId = plt.brep; end
            if nargin < 3, plotLabels = 1; end
            if nargin < 4, title = 'Geometry'; end

            % Retrieve geometry data
            pdeGeom = obj.m_pdeGeom;
            segmentMapping = obj.m_brepSegsMapping;

            % Plot geometry
            fig = figure(figId); clf(fig,'reset');
            pdegplot(pdeGeom);
            set(gcf, 'Name', title);
            xlabel('$x$', 'Interpreter', 'latex');
            ylabel('$y$', 'Interpreter', 'latex');

            hold on;

            % Unique segments
            uniqueSegments = unique(segmentMapping);

            if plotLabels == 1
                for i = 1:length(uniqueSegments)
                    brepSegment = uniqueSegments(i);
                    segIndices = find(segmentMapping == brepSegment);

                    midpoints = zeros(2, numel(segIndices));

                    for k = 1:numel(segIndices)
                        seg = segIndices(k);
                        type = pdeGeom(1,seg);

                        if type == 1  % Arc
                            xs = pdeGeom(2,seg); xe = pdeGeom(3,seg);
                            ys = pdeGeom(4,seg); ye = pdeGeom(5,seg);
                            xc = pdeGeom(8,seg); yc = pdeGeom(9,seg);

                            thetaStart = atan2(ys - yc, xs - xc);
                            thetaEnd = atan2(ye - yc, xe - xc);

                            angleDiff = mod(thetaEnd - thetaStart, 2*pi);
                            if angleDiff > pi
                                angleMid = thetaStart - angleDiff/2;
                            else
                                angleMid = thetaStart + angleDiff/2;
                            end

                            xmid = xc + pdeGeom(10,seg) * cos(angleMid);
                            ymid = yc + pdeGeom(10,seg) * sin(angleMid);

                        elseif type == 2  % Line
                            xmid = mean(pdeGeom(2:3,seg));
                            ymid = mean(pdeGeom(4:5,seg));
                        end

                        midpoints(:,k) = [xmid; ymid];
                    end

                    % Choose the segment closest to the centroid as label location
                    centroid = mean(midpoints,2);
                    distances = vecnorm(midpoints - centroid,2,1);
                    [~, idx] = min(distances);

                    xLabel = midpoints(1,idx);
                    yLabel = midpoints(2,idx);

                    text(xLabel, yLabel, num2str(brepSegment), 'Color', 'r', ...
                        'FontWeight', 'bold', 'HorizontalAlignment', 'center');
                end

                % Plot vertices
                p = obj.m_brep.vertices;
                plot(p(1,:),p(2,:),'b*');
                for i = 1:size(p,2)
                    text(p(1,i),p(2,i),num2str(i),'Color','b','VerticalAlignment','bottom');
                end
            end

            hold off;
            drawnow;
        end

        function plotPolygon(obj)
            % PLOTPOLYGON Function to plot the polygon
            % Input Arguments:
            %     obj - object containing the polygon data
            % Output Arguments:
            %     None, but creates a figure with the polygon plotted

            plt = PlotId;
            % Create a new figure for plotting
            fig = figure(plt.polygon);clf(fig,'reset');
            % Plot the geometry
            plot(polyshape(obj.m_polygon(1,:),obj.m_polygon(2,:)))
            title = 'Polygon';
            set(gcf, 'Name', title);
            xlabel('$x$')
            ylabel('$y$')
            axis("tight")
        end

        function drawArrow(~,start,stop,clr,thickness)
            % DRAWARROW Function to draw arrows between start and stop points
            %
            % Input Arguments:
            %     start - Nx2 matrix of starting points
            %     stop - Nx2 matrix of ending points
            %     clr - color of the arrow (optional)
            %     thickness - thickness of the arrow line (optional)

            % Set default color and thickness if not provided
            if (nargin == 3)
                clr = 'k';
                thickness = 1.0;
            end
            % Loop through each pair of start and stop points
            for i = 1:size(start,1)
                x0 = start(i,1);    y0 =start(i,2);
                x1 = stop(i,1);     y1 = stop(i,2);
                % Plot the line segment representing the arrow
                h = plot([x0 x1],[y0, y1],'Color',clr);
                set(h,'linewidth',thickness);
                hold on;
                % Calculate the direction vector for the arrow
                p = stop(i,:)-start(i,:);
                alpha = 0.4;  % Size of arrow head relative to the length of the vector
                beta = 0.35;  % Width of the base of the arrow head relative to the length
                % Calculate the coordinates for the arrow head
                hu = [x1-alpha*(p(1)+beta*(p(2)+eps)); x1; x1-alpha*(p(1)-beta*(p(2)+eps))];
                hv = [y1-alpha*(p(2)-beta*(p(1)+eps)); y1; y1-alpha*(p(2)+beta*(p(1)+eps))];
                % Plot the arrow head
                h = plot(hu(:),hv(:),'Color',clr);
                set(h,'linewidth',thickness);
            end
        end

        function pt = findPointinEmptyRegion(obj,x_limits,y_limits,dataPoints)
            % FINDPOINTINEMPTYREGION Function to find a point in an empty region
            %
            % Input Arguments:
            %     obj - object containing geometric data
            %     x_limits - limits for x-axis
            %     y_limits - limits for y-axis
            %     dataPoints - existing data points to avoid

            x_poly = [];
            y_poly = [];
            numPoints = 10; % Number of points to sample on the edge
            segments = 1:size(obj.m_pdeGeom,2);
            % Loop through each segment of the geometry
            for seg = segments
                % Sample points along the edges of the geometry
                x = linspace(min(obj.m_pdeGeom(2:3,seg)),max(obj.m_pdeGeom(2:3,seg)),numPoints);
                y = linspace(min(obj.m_pdeGeom(4:5,seg)),max(obj.m_pdeGeom(4:5,seg)),numPoints);
                x_poly = [x_poly x]; %#ok
                y_poly = [y_poly y]; %#ok
            end
            % Define a grid over the plot
            num_grid_points = 100; % Adjust for finer resolution
            x_grid = linspace(x_limits(1), x_limits(2), num_grid_points);
            y_grid = linspace(y_limits(1), y_limits(2), num_grid_points);
            [X_grid, Y_grid] = meshgrid(x_grid, y_grid);
            gridPts_x = X_grid(:);  gridPts_y = Y_grid(:);
            gridPts = [gridPts_x';gridPts_y'];
            % Check which grid points are inside the boundary representation
            in = obj.inBrep(gridPts);
            % Get the plotted data points
            x_data = [ dataPoints(1,:),x_poly];
            y_data = [ dataPoints(2,:),y_poly];
            % Compute the distance from each grid point to the closest data point
            min_distances = zeros(size(X_grid));
            for i = 1:numel(in)
                if ~in(i), continue;end
                % For each grid point, find the distance to all data points and take the minimum
                distances = sqrt((x_data - gridPts(1,i)).^2 + (y_data - gridPts(2,i)).^2);
                min_distances(i) = min(distances);
            end
            % Find the grid point with the maximum distance (i.e., furthest from existing data)
            [~, idx] = max(min_distances(:));
            % Get the corresponding coordinates of the empty space
            pt = [X_grid(idx);Y_grid(idx)];
        end

        %% Get methods
        function brep = getBrep(obj)
            % Method to retrieve the boundary representation (BRep) of the object
            brep = obj.m_brep;
        end

        function v = getBrepVertices(obj)
            % Method to retrieve the vertices of boundary representation (BRep) of the object
            v = obj.m_brep.vertices;
        end

        function v = getBrepSegments(obj)
            % Method to retrieve the segments of boundary representation (BRep) of the object
            v = obj.m_brep.vertices;
        end

        function numVertices = getNumberofBrepVertices(obj)
            % Method to retrieve the number of vertices of boundary representation (BRep) of the object
            numVertices = size(obj.m_brep.vertices,2);
        end

        function numSegments = getNumberofBrepSegments(obj)
            % Method to retrieve the number of vertices of boundary representation (BRep) of the object
            numSegments = size(obj.m_brep.segments,2);
        end

        function area = getArea(obj)
            % Method to retrieve the the precomputed area of the object
            area = obj.m_area;
        end

        function perimeter = getPerimeter(obj)
            % Method to retrieve the the precomputed perimeter of the object
            perimeter = obj.m_perimeter;
        end

        function pdeGeom = getPdeGeom(obj)
            % Method to retrieve the PDE geometry of the object
            pdeGeom = obj.m_pdeGeom;
        end

        function polygon = getPolygon(obj)
            % Method to retrieve the polygon representation of the object
            polygon = obj.m_polygon;
        end

        function numSegsPerArc = getNumSegsPerArc(obj)
            % Method to retrieve the number of segments per arc
            numSegsPerArc = obj.m_numSegsPerArc;
        end

        function brepSegsMapping = getBrepSegsMapping(obj)
            % Method to retrieve the mapping of BRep segments
            brepSegsMapping = obj.m_brepSegsMapping;
        end

        function numBndrySegs = getNumBndrySegs(obj)
            % Method to retrieve the number of boundary segments
            numBndrySegs = obj.m_numBndrySegs;
        end

        function segLengths = getSegLengths(obj)
            % Method to retrieve the lengths of segments
            segLengths = obj.m_segLengths;
        end

        function modelScale = getModelScale(obj)
            % Method to retrieve the scale of the model
            modelScale = obj.m_modelScale;
        end

        function boundingBox = getBoundingBox(obj)
            % Method to retrieve the bounding box of the object
            boundingBox = obj.m_boundingBox;
        end

        function boxSizes = getBoxSizes(obj)
            % Method to retrieve the sizes of the boxes
            boxSizes = obj.m_boxSizes;
        end

        function center = getCenter(obj)
            % Method to retrieve the center of the object
            center = obj.m_center;
        end

        function adjustFigScale(~)
            axis tight;
            v = axis;
            xRange = v(2)-v(1);
            yRange = v(4)-v(3);
            deltaX = 0.1*xRange;
            deltaY = 0.1*yRange;
            v = [v(1)-deltaX v(2)+deltaX  v(3)-deltaY v(4)+deltaY];
            axis(v);
        end
    end
end