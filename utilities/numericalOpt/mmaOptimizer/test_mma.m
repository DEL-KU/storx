n = 3;  m = 2;
x = [4;3;2];
xold1 = x;
xold2   = x;
xmin = zeros(n,1);  xmax = 5*ones(n,1);
mma = mmaOptimizer(n,m,0,1000,1,xold);

low     = xmin;
upp     = xmax;
c       = [1000  1000]';
d       = [1  1]';
a      = 1;
ai       = [0  0]';

change = 1;
for iter = 1:20
    if change < 0.001,break;end
    f0val = x(1)^2 + x(2)^2 + x(3)^2;
%
df0dx = [2*x(1)
	 2*x(2)
	 2*x(3)];
%
fval  = [(x(1)-5)^2+(x(2)-2)^2+(x(3)-1)^2-9
	 (x(1)-3)^2+(x(2)-4)^2+(x(3)-3)^2-9];
%
dfdx  = 2*[x(1)-5  x(2)-2  x(3)-1
	   x(1)-3  x(2)-4  x(3)-3];


% [x,~,~,~,~,~,~,~,~,low,upp] = ...
%                 mmasub(m,n,iter,x,xmin,xmax,xold1,xold2, ...
%                 f0val,df0dx,fval,dfdx,low,upp,a,ai,c,d);
    % 
    [mma,x] = mma.update(x,f0val,df0dx,fval,dfdx,xmin,xmax);
    change = max(abs(x - xold1));
    xold1 = x;
xold2 = xold1;
    disp(['iter: ' num2str(iter) ', change:' num2str(change) ', f(x): ', num2str(f0val) ', x: ', num2str(x')])
end



