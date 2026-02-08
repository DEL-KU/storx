%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for exporting .DXF file from a design pseudo-density matrix          %
% Modified "export.m" code in:                                              %
% NAVIER-STOKES TOPOLOGY OPTIMISATION CODE, MAY 2022                        %
% COPYRIGHT (c) 2022, J ALEXANDERSEN. BSD 3-CLAUSE LICENSE                  %
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


function dxf_filename = export_dxf(field, level, bbox, hx, hy, example_name)
% exportDXF_fromContourcOnBox
% field   : ny-by-nx scalar field (e.g., -xPhys)
% level   : contour level (e.g., -0.5)
% bbox    : [Xmin Xmax; Ymin Ymax]
% hx, hy  : grid spacing in x/y
% filename: output .dxf

arguments
    field double
    level (1,1) double
    bbox (2,2) double
    hx (1,1) double
    hy (1,1) double
    example_name (1,:) char
end

Xmin = bbox(1,1) - hx;
Xmax = bbox(1,2) + hx;
Ymin = bbox(2,1) - hy; 
Ymax = bbox(2,2) + hy;

field_padded = zeros(size(field)+2);
field_padded(2:end-1,2:end-1) = field;

C = contourc(field_padded,[level level]);
C = unique(C','rows','stable')';

C = C-2.5; C(1,:) = C(1,:)*hx; C(2,:) = C(2,:)*hy;

h = hypot(hx,hy);

ind1 = intersect(find(C(1,:)<Xmax+hx),find(C(2,:)<Ymax+hy));
ind2 = intersect(find(C(1,:)>Xmin-hx),find(C(2,:)>Ymin-hy));
ind = intersect(ind1,ind2);
X = []; Y = []; c = 0;
indused = [];
for i = 1:(length(ind))
    x0 = C(1,ind(i)); y0 = C(2,ind(i));
    redind = ind([1:max(1,i-5) i+1:min(length(ind),i+2)]);
    redind = setdiff(redind,indused(max(1,c-10):end));
    R = sqrt( (C(1,redind)-x0).^2 + (C(2,redind)-y0).^2 );
    [minR,idx] = min(R);
    if (minR < 2*h && minR > 0)
        c = c + 1;
        indused(c) = i; %#ok
        X(c,:) = [x0 C(1,redind(idx))]; %#ok
        Y(c,:) = [y0 C(2,redind(idx))]; %#ok
    end
end
% Ensure .dxf extension
[~,~,ext] = fileparts(example_name);
if isempty(ext)
    dxf_filename = [example_name '.dxf'];
elseif strcmpi(ext,'.dxf')
    dxf_filename = filename;
else
    dxf_filename = [filename '.dxf']; % if user gave some other extension, append .dxf
end

fid = fopen(dxf_filename,'w');

fprintf(fid,'0\nSECTION\n2\nENTITIES\n0\n');

% ---- close loops by chaining segments ----
nseg = size(X,1);
used = false(nseg,1);

x1 = X(:,1); 
y1 = Y(:,1);
x2 = X(:,2);
y2 = Y(:,2);

% choose a tolerance based on your grid (adjust if needed)
rx = max([x1; x2]) - min([x1; x2]);
ry = max([y1; y2]) - min([y1; y2]);
tol = 1e-6 * max(rx, ry);
if tol == 0, tol = 1e-9; end

while any(~used)

    % start a new chain from first unused segment
    i0 = find(~used, 1, 'first');
    used(i0) = true;

    loop = [x1(i0) y1(i0); x2(i0) y2(i0)];  % ordered vertices
    cur  = loop(end,:);

    % grow chain by matching endpoints
    while true
        ds = hypot(x1 - cur(1), y1 - cur(2));  % match next start
        de = hypot(x2 - cur(1), y2 - cur(2));  % match next end (needs flip)

        ds(used) = inf;
        de(used) = inf;

        [mS, iS] = min(ds);
        [mE, iE] = min(de);

        if mS <= tol
            % connect to segment as-is: (x1,y1)->(x2,y2)
            used(iS) = true;
            nxt = [x2(iS) y2(iS)];
        elseif mE <= tol
            % connect by flipping segment: (x2,y2)->(x1,y1)
            used(iE) = true;
            nxt = [x1(iE) y1(iE)];
        else
            % no continuation found
            break;
        end

        loop(end+1,:) = nxt; %#ok<AGROW>
        cur = nxt;

        % if we returned to the start, stop
        if hypot(cur(1) - loop(1,1), cur(2) - loop(1,2)) <= tol
            break;
        end
    end

    % force closure if needed
    if hypot(loop(end,1) - loop(1,1), loop(end,2) - loop(1,2)) > tol
        loop(end+1,:) = loop(1,:); %#ok  % add closing vertex
    end

    % write this loop as LINE entities
    for j = 1:(size(loop,1)-1)
        fprintf(fid,'LINE\n8\n0\n');
        fprintf(fid,'10\n%.6f\n20\n%.6f\n30\n%.6f\n', loop(j,1),   loop(j,2),   0);
        fprintf(fid,'11\n%.6f\n21\n%.6f\n31\n%.6f\n', loop(j+1,1), loop(j+1,2), 0);
        fprintf(fid,'0\n');
    end
end

fprintf(fid,'ENDSEC\n0\nEOF\n');
fclose(fid);


end
