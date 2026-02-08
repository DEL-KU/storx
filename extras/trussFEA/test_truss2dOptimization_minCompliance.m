% truss optimize examples
clc; close all; format compact; format short

problem = 1;
trussClass = 'truss2dMinComplianceVolumeConstraint';

if (problem == 0) % two bar truss problem with analytical solution
    theta = pi/4;% bar orientation
    beta = pi/3;%angle of force
    P = 10;
    xy = [0 -cos(theta) cos(theta); 0 sin(theta) sin(theta)];
    connectivity = [1 2; 1 3]';
    t = feval(trussClass,xy,connectivity);%#ok<FVAL> % initialize model
    t = t.assignE(2e11); 
    t = t.assignA(1e-6);
    t = t.fixXofNodes([2 3]);
    t = t.fixYofNodes([2 3]);
    t = t.applyForce(1,[P*cos(beta); -P*sin(beta)]);
elseif (problem == 1)
    xy = [0.5 1.5 1.0 0 2.0; 1 1 0 0 0];% nodes
    connectivity = [1 2; 1 3; 1 4; 2 3; 2 5; 3 4; 3 5]'; % connectivity  
    t = feval(trussClass,xy,connectivity);%#ok<FVAL> % initialize model
    t = t.assignE(2e11);
    t = t.assignA(1e-6); % for all members 
    t = t.assignYieldStress(100e6); 
    t = t.fixXofNodes([4 5]);
    t = t.fixYofNodes([4 5]);
    t = t.applyForce(1,[1;-2]);
    t = t.applyForce(2,[2;0]);
elseif (problem == 2)
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
    t = feval(trussClass,xy,connectivity);% initialize model
    t = t.assignE(2e11);
    t = t.assignA(1e-6); % for all members
    
    t = t.fixXofNodes(1:M);
    t = t.fixYofNodes(1:M);
    t = t.applyForce(N*M-round((M-1)),[0;-1]);
end
t.plot(1,1);
tic
t = t.optimize();
toc
optimizedArea = t.myArea;
t.plotDeformed(3);
disp('Initial Volume: '); disp(t.myInitialVolume)
disp('Final Volume: '); disp(t.myFinalVolume);
disp('Initial Compliance: '); disp(t.myInitialCompliance)
disp('Final Compliance: '); disp(t.myFinalCompliance);
disp('Initial Stress/Yield Strength: '); disp(max(abs(t.myInitialStress/t.myYieldStress)))
disp('Final Stress/Yield Strength: '); disp(max(abs(t.myFinalStress/t.myYieldStress)));