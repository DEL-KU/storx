function [dl2, arcReport] = pdegeom_fixArcsProjectToCircle(dl, tol)
% For each arc segment (type=1), enforce endpoints lie on circle (xc,yc,R)
% by projecting (x1,y1) and (x2,y2) to the circle.
%
% Also optionally recomputes R from average of endpoint radii if inconsistent.

    arguments
        dl double
        tol double = NaN
    end

    dl2 = dl;

    typ = dl(1,:);
    x1  = dl(2,:); x2 = dl(3,:);
    y1  = dl(4,:); y2 = dl(5,:);

    PALL = [x1(:),y1(:); x2(:),y2(:)];
    if isnan(tol)
        Lx = max(PALL(:,1)) - min(PALL(:,1));
        Ly = max(PALL(:,2)) - min(PALL(:,2));
        scl = max([Lx,Ly,1]);
        tol = 1e-10*scl;
    end

    arcIds = find(typ==1);
    arcReport = struct('seg',{},'d1',{},'d2',{},'Rold',{},'Rnew',{});

    if isempty(arcIds), return; end
    if size(dl,1) < 10
        error('Arc segments present but dl has <10 rows (need center+radius).');
    end

    for i = arcIds(:)'
        xc = dl2(8,i); yc = dl2(9,i); R = dl2(10,i);

        d1 = hypot(dl2(2,i)-xc, dl2(4,i)-yc);
        d2 = hypot(dl2(3,i)-xc, dl2(5,i)-yc);

        % If R is inconsistent with endpoints, choose a better R
        Rnew = R;
        if abs(d1-R) > 1e2*tol || abs(d2-R) > 1e2*tol
            Rnew = 0.5*(d1+d2);  % average radius from both endpoints
            dl2(10,i) = Rnew;
        end

        % Project endpoints onto circle of radius Rnew
        v1 = [dl2(2,i)-xc, dl2(4,i)-yc];
        v2 = [dl2(3,i)-xc, dl2(5,i)-yc];

        if norm(v1) > 0
            v1 = (Rnew/norm(v1))*v1;
            dl2(2,i) = xc + v1(1);
            dl2(4,i) = yc + v1(2);
        end
        if norm(v2) > 0
            v2 = (Rnew/norm(v2))*v2;
            dl2(3,i) = xc + v2(1);
            dl2(5,i) = yc + v2(2);
        end

        arcReport(end+1) = struct('seg',i,'d1',d1,'d2',d2,'Rold',R,'Rnew',Rnew); %#ok<AGROW>
    end
end
