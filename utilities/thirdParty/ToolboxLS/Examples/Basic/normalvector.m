function [ normalx, normaly ] = normalvector(grid, data)

% Qian Ye, 2/15/2019
%---------------------------------------------------------------------------
% Get the first and second derivative terms.
[ second, first ] = hessianSecond(grid, data);

%---------------------------------------------------------------------------
% Compute gradient magnitude.
gradMag2 = first{1}.^2;
for i = 2 : grid.dim
  gradMag2 = gradMag2 + first{i}.^2;
end
gradMag = sqrt(gradMag2);
normalx = zeros(size(data));
normaly =zeros(size(data));
normalx= first{1};
normaly=first{2};
normalx= first{1}./(gradMag+eps);
normaly=first{2}./(gradMag+eps);
