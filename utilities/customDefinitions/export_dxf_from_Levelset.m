function export_dxf_from_Levelset(phiF, xF, yF, level, filename, minPts)
% Export DXF polylines from iso-contour of phiF at "level"
% Solid is assumed phi <= level (your contourf([-Inf 0]) convention).
% Padding with "void" makes contours close even when solid touches the bbox.

arguments
    phiF double
    xF double
    yF double
    level (1,1) double
    filename (1,:) char
    minPts (1,1) double {mustBeNonnegative} = 0
end

% Ensure vectors
xF = xF(:)';     % 1-by-Nx
yF = yF(:)';     % 1-by-Ny

% Grid step (assumes uniform-ish)
dx = mean(diff(xF));
dy = mean(diff(yF));

% Choose a padding value that is definitely "void" for phi<=level solids
v = phiF(isfinite(phiF));
if isempty(v), error('phiF has no finite values.'); end
padVal = max(v) + abs(max(v)-min(v)) + 1;   % safely > level in practice

% Treat NaNs/Inf as void too (same intent as your old inDom mask)
phi = phiF;
phi(~isfinite(phi)) = padVal;

% Pad phi with a 1-cell void border
phiP = padarray(phi, [1 1], padVal, 'both');

% Pad coordinate vectors consistently
xP = [xF(1)-dx, xF, xF(end)+dx];
yP = [yF(1)-dy, yF, yF(end)+dy];

% Get contour matrix from contourc on padded grid
C = contourc(xP, yP, phiP, [level level]);

% Output name
[~,~,ext] = fileparts(filename);
if isempty(ext), outname = [filename '.dxf'];
elseif strcmpi(ext,'.dxf'), outname = filename;
else, outname = [filename '.dxf'];
end

fid = fopen(outname,'w');  assert(fid>0, 'Cannot write %s', outname);
fprintf(fid,'0\nSECTION\n2\nENTITIES\n');

% Parse contourc: [level; npts] header then npts columns of [x;y]
j = 1;
while j < size(C,2)
    npts = C(2,j);
    P = C(:, j+1 : j+npts);
    j = j + npts + 1;

    if minPts > 0 && npts < minPts
        continue
    end

    x = P(1,:).';
    y = P(2,:).';

    % Ensure closed
    if x(1)~=x(end) || y(1)~=y(end)
        x(end+1,1) = x(1);
        y(end+1,1) = y(1);
    end

    n = numel(x)-1;
    fprintf(fid,'0\nLWPOLYLINE\n8\nCONTOUR\n');
    fprintf(fid,'90\n%d\n', n);
    fprintf(fid,'70\n1\n'); % closed
    for p = 1:n
        fprintf(fid,'10\n%.6f\n20\n%.6f\n', x(p), y(p));
    end
end

fprintf(fid,'0\nENDSEC\n0\nEOF\n');
fclose(fid);
end
