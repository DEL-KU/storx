clc; close all;
clear; figure
format compact;
format long
clear t
problem = 7;
if (problem == 1) % verify against force balance
    xy = [0 0.5 0.25; 0 0 0.5];% nodes
    connectivity = [1 3; 2 3]'; % connectivity
    t = truss2d(xy,connectivity);
    t = t.assignE(2e11); % for all members  
    t = t.assignA(1e-6); % for all members 
    t = t.fixXofNodes([1 2]);
    t = t.fixYofNodes([1 2]);
    t = t.applyForce(3,[100; 25]);
elseif (problem == 2)
    xy = [0  1.25 2 1.25; 0 0 0 1];% nodes
    connectivity = [1 2; 1 4; 2 3; 2 4; 3 4]'; % connectivity
    t = truss2d(xy,connectivity);
    t = t.assignE(2e11); % for all bars  
    t = t.assignA(1e-6);
    t = t.assignA(2*1e-6,4);
    t = t.fixXofNodes(1);
    t = t.fixYofNodes([1 3]);
    t = t.applyForce(4,[0; -100]);
elseif (problem == 3)
    xy = [0 0 -1/sqrt(2) 1/sqrt(2);
          0 1 -1/sqrt(2) -1/sqrt(2) ];% nodes
    connectivity = [1 2; 1 3;  1 4]'; % connectivity
    t = truss2d(xy,connectivity);
    t = t.assignE(2e11); % for all members
     %t = t.assignE(0.01,2); 
    t = t.assignA(1e-6);
    t = t.fixXofNodes([2 3 4]);
    t = t.fixYofNodes([2 3 4]);
    t = t.applyForce(1,[4; 2]);
elseif (problem == 4)
    xy = [0.5 1.5 1.0 0 2.0; 1 1 0 0 0];% nodes
    connectivity = [1 2; 1 3; 1 4; 2 3; 2 5; 3 4; 3 5]'; % connectivity  
    t = truss2d(xy,connectivity);
    t = t.assignE(2e11);
    t = t.assignA(1e-6); % for all members 
    t = t.fixXofNodes([4 5]);
    t = t.fixYofNodes([4 5]);
    t = t.applyForce(1,[1;-2]);
    t = t.applyForce(2,[2;0]);
elseif (problem == 5)
    xy = [0.5 1.5 1.0 0 2.0; 1 1 0 0 0];% nodes
    connectivity = [1 2; 1 3; 1 4; 2 3; 2 5; 3 4; 3 5]'; % connectivity  
    t = truss2d(xy,connectivity);
    t = t.assignE(2e11);
    t = t.assignA(1e-6); % for all members 
    t = t.fixXofNodes([4 5]);
    t = t.fixYofNodes([4 5]);
    t = t.applyForce(3,[0;-1]);
elseif (problem == 6)
    xy = [-1 1 -2 2 2 -2; 0 0 -1 -1 1 1];% nodes
    connectivity = [1 2; 1 3; 1 6; 2 4; 2 5]'; % connectivity
    t = truss2d(xy,connectivity);
    t = t.assignE(1000*2,1);
    t = t.assignE(1000*sqrt(2),[2 3 4 5]);
    t = t.assignA(1); % for all members 
    t = t.fixXofNodes([3 4 5 6]);
    t = t.fixYofNodes([3 4 5 6]);
    t = t.applyForce(1,[1;-2]);
    t = t.applyForce(1,[100; -25]);
elseif (problem == 7) % box truss
    xy = [0  1 1 0; 0 0 1 1];% nodes
    connectivity = [1 2; 2 3; 3 4; 4 1; 1 3 ]'; % connectivity
    t = truss2d(xy,connectivity);
    t = t.assignE(2e11); % for all members  
    t = t.assignA(1e-6); % for all members 
    t = t.fixXofNodes([1 2]);
    t = t.fixYofNodes([1 2]);
    t = t.applyForce(3,[100; 25]);
elseif (problem == 8)
    % grid of trusses
    M = 6;
    N = 2*M-1;
    LengthScale = 1e-6;
    [X,Y] = meshgrid(1:1:N, 1:1:M);
    xy = zeros(2,N*M);
    xy(1,:) = LengthScale*X(:);
    xy(2,:) = LengthScale*Y(:);
    connectivity = zeros(2,(N-1)*M + (M-1)*N + 2*(N-1)*(M-1));
    count = 1;
    % Vertical bars
    for i = 1:N
        for j = 1:M-1
            connectivity(1,count) =  (i-1)*(M) + j;
            connectivity(2,count) =  (i-1)*(M) + j+1;
            count  = count+1;
        end
    end
    % Horizontal bars
    for j = 1:M
        for i = 1:N-1
            connectivity(1,count) =  (i-1)*(M) + j;
            connectivity(2,count) =  (i)*(M) + j;
            count  = count+1;
        end
    end
    % Cross bars
    for i = 1:N-1
        for j = 1:M-1
            connectivity(1,count) =  (i-1)*(M) + j;
            connectivity(2,count) =  (i)*(M) + j+1;
            count  = count+1;
            connectivity(1,count) =  (i-1)*(M) + j+1;
            connectivity(2,count) =  (i)*(M) + j;
            count  = count+1;
        end
    end
    t = truss2d(xy,connectivity);% initialize model
    t = t.assignE(2e11);
    t = t.assignA(1e-6); % for all members
    t = t.fixXofNodes(1:M);
    t = t.fixYofNodes(1:M);
    t = t.applyForce(N*M-round((M-1)),[0;-1]);
elseif (problem == 9) % hanging node
    xy = [0 0.5 0.25 0.5; 0 0 0.5 0.5];% nodes
    connectivity = [1 3; 2 3]'; % connectivity
    t = truss2d(xy,connectivity);
    t = t.assignE(2e11); % for all members  
    t = t.assignA(1e-6); % for all members 
    t = t.fixXofNodes([1 2]);
    t = t.fixYofNodes([1 2]);
    t = t.applyForce(3,[100; 25]);  
elseif (problem == 10) % hanging bar
    xy = [0 0.5 0.25 0.5; 0 0 0.5 0.5];% nodes
    connectivity = [1 3; 2 3; 3 4]'; % connectivity
    t = truss2d(xy,connectivity);
    t = t.assignE(2e11); % for all members  
    t = t.assignA(1e-6); % for all members 
    t = t.fixXofNodes([1 2]);
    t = t.fixYofNodes([1 2]);
    t = t.applyForce(3,[100; 25]);  
elseif (problem == 11) % another example of hanging bar
    xy = [0.5 1.5 1.0 0 2.0 1; 1 1 0 0 0 0.5];% nodes
    connectivity = [1 2; 1 3; 1 4; 2 3; 2 5; 3 4; 3 5; 1 6]'; % connectivity  
    t = truss2d(xy,connectivity);
    t = t.assignE(2e11);
    t = t.assignA(1e-6); % for all members 
    t = t.fixXofNodes([4 5]);
    t = t.fixYofNodes([4 5]);
    t = t.applyForce(1,[1;-2]);
    t = t.applyForce(2,[2;0]);
elseif (problem == 12) % insufficient node connectivity
    xy = [0 0.5 0.25 0.375; 0 0 0.5 0.25];% nodes
    connectivity = [1 3; 3 4; 4 2]'; % connectivity
    t = truss2d(xy,connectivity);
    t = t.assignE(2e11); % for all members  
    t = t.assignA(1e-6); % for all members 
    t = t.fixXofNodes([1 2]);
    t = t.fixYofNodes([1 2]);
    t = t.applyForce(3,[100; 25]);  
elseif (problem == 13) % insufficient node connectivity
    xy = [0 0.5 0.25 0.4; 0 0 0.5 0.3];% nodes
    connectivity = [1 3; 3 4; 4 2]'; % connectivity
    t = truss2d(xy,connectivity);
    t = t.assignE(2e11); % for all members  
    t = t.assignA(1e-6); % for all members 
    t = t.fixXofNodes([1 2]);
    t = t.fixYofNodes([1 2]);
    t = t.applyForce(3,[100; 25]);

end
t.plot(1,1);
t = t.assemble();
t = t.solve();
t.plotDeformed(2);
t = t.computeStresses();
disp('Displacement:')
t.myUV

disp('Stress:')
t.myStress