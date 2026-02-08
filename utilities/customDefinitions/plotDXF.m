%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code for plotting .DXF file                                               %
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

function S = plotDXF(dxfFile, varargin)
% plotDXF  Plot polylines/lines from a DXF file.
% Supports: LWPOLYLINE, POLYLINE/VERTEX/SEQEND, LINE
%
% Usage:
%   plotDXF("shape.dxf");
%   S = plotDXF("shape.dxf","Fill",true);

p = inputParser;
addParameter(p,'Ax',[],@(h) isempty(h) || isgraphics(h,'axes'));
addParameter(p,'Fill',false,@(x) islogical(x) || isnumeric(x));
addParameter(p,'ShowPoints',false,@(x) islogical(x) || isnumeric(x));
parse(p,varargin{:});
ax = p.Results.Ax;
doFill = logical(p.Results.Fill);
showPts = logical(p.Results.ShowPoints);

if isempty(ax)
    figure('Color','w'); ax = axes; hold(ax,'on');
else
    hold(ax,'on');
end
axis(ax,'equal'); grid(ax,'on');

txt = readlines(dxfFile);
txt = strip(txt);
txt(txt=="") = [];

i = 1;
polys = {};
polyLayers = {};
lines = {};
lineLayers = {};

while i <= numel(txt)-1
    code = str2double(txt(i));
    val  = txt(i+1);

    if ~isnan(code) && code==0
        ent = char(val);

        switch upper(ent)
            case 'LWPOLYLINE'
                % parse until next "0 <ENTITY>"
                i = i + 2;
                layer = "0";
                nVert = NaN;
                closed = false;
                X = []; Y = [];

                while i <= numel(txt)-1
                    c = str2double(txt(i));
                    v = txt(i+1);

                    if ~isnan(c) && c==0
                        break; % next entity
                    end

                    switch c
                        case 8
                            layer = string(v);
                        case 90
                            nVert = str2double(v); %#ok<NASGU>
                        case 70
                            flags = str2double(v);
                            closed = bitand(flags,1)~=0;
                        case 10
                            X(end+1,1) = str2double(v); %#ok<AGROW>
                        case 20
                            Y(end+1,1) = str2double(v); %#ok<AGROW>
                        otherwise
                            % ignore other codes (bulge 42, widths 40/41, etc.)
                    end

                    i = i + 2;
                end

                if ~isempty(X) && numel(X)==numel(Y)
                    P = [X(:) Y(:)];
                    if closed && ~isequal(P(1,:),P(end,:))
                        P(end+1,:) = P(1,:);
                    end
                    polys{end+1} = P; %#ok<AGROW>
                    polyLayers{end+1} = layer; %#ok<AGROW>
                end

            case 'POLYLINE'
                % classic POLYLINE with VERTEX records until SEQEND
                i = i + 2;
                layer = "0";
                closed = false;

                % read POLYLINE header fields
                while i <= numel(txt)-1
                    c = str2double(txt(i));
                    v = txt(i+1);
                    if ~isnan(c) && c==0
                        break;
                    end
                    if c==8,  layer = string(v); end
                    if c==70
                        flags = str2double(v);
                        closed = bitand(flags,1)~=0;
                    end
                    i = i + 2;
                end

                % now expect many VERTEX blocks
                X = []; Y = [];
                while i <= numel(txt)-1
                    c0 = str2double(txt(i));
                    v0 = txt(i+1);
                    if ~isnan(c0) && c0==0 && strcmpi(v0,'VERTEX')
                        i = i + 2;
                        vx = NaN; vy = NaN;
                        while i <= numel(txt)-1
                            c = str2double(txt(i));
                            v = txt(i+1);
                            if ~isnan(c) && c==0
                                break;
                            end
                            if c==10, vx = str2double(v); end
                            if c==20, vy = str2double(v); end
                            i = i + 2;
                        end
                        if ~isnan(vx) && ~isnan(vy)
                            X(end+1,1)=vx; Y(end+1,1)=vy; %#ok<AGROW>
                        end
                    elseif ~isnan(c0) && c0==0 && strcmpi(v0,'SEQEND')
                        i = i + 2;
                        break;
                    else
                        % unknown record inside polyline; advance
                        i = i + 2;
                    end
                end

                if ~isempty(X) && numel(X)==numel(Y)
                    P = [X(:) Y(:)];
                    if closed && ~isequal(P(1,:),P(end,:))
                        P(end+1,:) = P(1,:);
                    end
                    polys{end+1} = P; %#ok<AGROW>
                    polyLayers{end+1} = layer; %#ok<AGROW>
                end

            case 'LINE'
                i = i + 2;
                layer = "0";
                x1=NaN;y1=NaN;x2=NaN;y2=NaN;
                while i <= numel(txt)-1
                    c = str2double(txt(i));
                    v = txt(i+1);
                    if ~isnan(c) && c==0
                        break;
                    end
                    switch c
                        case 8,  layer = string(v);
                        case 10, x1 = str2double(v);
                        case 20, y1 = str2double(v);
                        case 11, x2 = str2double(v);
                        case 21, y2 = str2double(v);
                    end
                    i = i + 2;
                end
                if all(~isnan([x1 y1 x2 y2]))
                    lines{end+1} = [x1 y1; x2 y2]; %#ok<AGROW>
                    lineLayers{end+1} = layer; %#ok<AGROW>
                end

            otherwise
                i = i + 2;
        end
    else
        i = i + 1;
    end
end

% Plot polylines
for k = 1:numel(polys)
    P = polys{k};
    plot(ax, P(:,1), P(:,2), '-', 'LineWidth', 1.0);
    if showPts
        plot(ax, P(:,1), P(:,2), '.', 'MarkerSize', 6);
    end
    if doFill && size(P,1) >= 3
        try
            pg = polyshape(P(:,1), P(:,2), 'Simplify', false);
            if pg.NumRegions > 0
                plot(ax, pg, 'FaceAlpha', 0.15, 'EdgeAlpha', 0.9);
            end
        catch
            % ignore non-simple polygons
        end
    end
end

% Plot lines
for k = 1:numel(lines)
    P = lines{k};
    plot(ax, P(:,1), P(:,2), '-', 'LineWidth', 1.0);
end

title(ax, sprintf('DXF: %s', dxfFile), 'Interpreter','none');
xlabel(ax,'X'); ylabel(ax,'Y');

S.polylines = polys;
S.polyLayers = polyLayers;
S.lines = lines;
S.lineLayers = lineLayers;
end
