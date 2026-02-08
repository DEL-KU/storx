function plot2DContour(field,value)
clf;
[nely,nelx] = size(field);
if (max(field) - min(field) < 1e-10)
    field = field + 1e-8*rand(size(field)); % when is constant
end
fill([1 nelx nelx 1],[1 1 nely nely],'b'); hold on;
contourf(-field,[-value -value]); axis('equal'); axis tight;axis off;pause(0.01);
nElemsAbove = sum(field(:)> value);
title(['Vol fraction: ' num2str(nElemsAbove/(nelx*nely))]);