function dl10 = pdegeom_normalizeTo10Rows(dl)
% Ensure decomposed geometry matrix is 10xN, with proper padding:
% Row 1: type (2 line, 1 circle arc)
% Row 2-5: x1 x2 y1 y2
% Row 6-7: left right
% Row 8-9: xc yc (arc only)
% Row 10: R (arc only)
% For lines: rows 8-10 set to 0.

    if size(dl,1) < 7
        error('dl must have at least 7 rows (needs left/right labels).');
    end

    nSeg = size(dl,2);
    dl10 = zeros(10,nSeg);

    % Copy common rows
    dl10(1,:) = dl(1,:);
    dl10(2:7,:) = dl(2:7,:);

    typ = dl10(1,:);

    % If arcs exist, require center+R in source dl
    if any(typ==1)
        if size(dl,1) < 10
            error('Arc segments present but dl has <10 rows (need center rows 8-9 and radius row 10).');
        end
        dl10(8,:)  = dl(8,:);
        dl10(9,:)  = dl(9,:);
        dl10(10,:) = dl(10,:);
    end

    % For lines, zero out arc rows explicitly
    isLine = (typ==2);
    dl10(8,isLine)  = 0;
    dl10(9,isLine)  = 0;
    dl10(10,isLine) = 0;
end
