clc; close all;
clear; figure
format compact;
format long
clear t

L = 1.0;
theta = pi/4;
beta = pi/3;
P = 1;
E = 200e9;
A1 = 1e-6 ;
A2 = 1e-6 ;
xy = [0 -L*cos(theta) L*cos(theta); ...
    0 L*sin(theta) L*sin(theta)];
connectivity = [1 2; 1 3]';

t = truss2d(xy,connectivity);
t = t.assignE(E);
t = t.assignA(A1,1);
t = t.assignA(A2,2);
t = t.fixXofNodes([2 3]);
t = t.fixYofNodes([2 3]);
t = t.applyForce(1,[P*cos(beta); -P*sin(beta)]);

%% Exact Solution
Fx = P*cos(beta)
Fy = -P*sin(beta)
k = E*A1/L;
disp(['Exact Solution:'])
uExact = Fx/(2*k*cos(theta)*cos(theta))
vExact = Fy/(2*k*sin(theta)*sin(theta))


%% Different Ai
% uExact = P*L/(2*E*cos(theta)*sin(2*theta))*(sin(beta+theta)/A1-sin(beta-theta)/A2);
% vExact = -P*L/(2*E*sin(theta)*sin(2*theta))*(sin(beta+theta)/A1+sin(beta-theta)/A2);
% disp([uExact vExact])

t.plot(1,1);
t = t.assemble();
t = t.solve();
t.plotDeformed(2);
t = t.computeStresses();

disp('Numerical Solution:')
t.myUV