% truss optimize examples
clc; close all; format compact; format short


problem =1;
nSamples = 1000;
if (problem == 1)
    xy = [0.5 1.5 1.0 0 2.0; 1 1 0 0 0];% nodes
    connectivity = [1 2; 1 3; 1 4; 2 3; 2 5; 3 4; 3 5]'; % connectivity
    t = truss2d(xy,connectivity);% initialize model
    t = t.assignE(2e11);
    t = t.assignA(1e-6); % for all members
    %t = t.assignYieldStress(100e6);
    t = t.fixXofNodes([4 5]);
    t = t.fixYofNodes([4 5]);
    t = t.applyForce(1,[1;-2]);
    t = t.applyForce(2,[2;0]);

    t = t.assemble();
    t = t.solve();

    t.plot(0);


    [relV,relJ] = trussSamplingRandom(t,nSamples);
end