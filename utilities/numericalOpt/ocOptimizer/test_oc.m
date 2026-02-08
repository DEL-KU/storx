close all; clear; clc;format long

[X,Y] = meshgrid(-2:0.1:2,-2:0.1:2);
a = 2;
b = 0.5;
fx = exp(a*(X-b))-a*(X-b) - 1 + exp(a*(Y-b))-a*(Y-b) - 1;

surf(X,Y,fx); alpha 0.5; hold on

contour(X,Y,fx,100)



%zlim([-2 5])

n = 2;
m = 1;
x = [-1.5;1.25];

xold = x;
xmin = 0;
xmax = 5;
gmax = -1;
oc =   ocOptimizer(0,1000,0.1,0,5,gmax);
disp('Test OC');
ch = 1;
for iter = 1:2000
    if ch < 2e-4,break;end
    fx = sum(x.^2);
    dfdx = [2*x(1)
        2*x(2)];
    
    
    [oc,x] = oc.update(x,-dfdx);
    
    ch = max(abs(x - xold));
    
    xold = x;
    
    fprintf('iter %d:\n',iter) ;
    fprintf('f =  %f\n', fx);
    fprintf('ch =  %f\n', ch);
    disp('g:');
    disp(sum(x))
    disp('x:')
    disp(x')
    
    
end