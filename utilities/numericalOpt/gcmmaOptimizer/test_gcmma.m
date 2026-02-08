n = 3;  m = 2;
x = [4;3;2];
xold = x;
xmin = zeros(n,1);  xmax = 5*ones(n,1);
gcmma = gcmmaOptimizer(n,m,0,1000,1,xold);
change = 1;
for iter = 1:50
    if change < 0.001,break;end
    fx = x(1)^2 + x(2)^2 + x(3)^2;
    dfdx = 2*x;
    gx = [(x(1)-5)^2 + (x(2)-2)^2 + (x(3)-1)^2 - 9
        (x(1)-3)^2 + (x(2)-4)^2 + (x(3)-3)^2 - 9];
    dgdx  = 2*[x(1)-5  x(2)-2  x(3)-1
        x(1)-3  x(2)-4  x(3)-3];

    [gcmma,xnew] = gcmma.outerUpdate(x,fx,dfdx,gx,dgdx,xmin,xmax);
    fxnew = xnew(1)^2 + xnew(2)^2 + xnew(3)^2;
    gxnew = [(xnew(1)-5)^2 + (xnew(2)-2)^2 + (xnew(3)-1)^2 - 9
        (xnew(1)-3)^2 + (xnew(2)-4)^2 + (xnew(3)-3)^2 - 9];

    conserv = gcmma.conCheck(fxnew,gxnew);
    for innerIter = 1:15
        if (conserv == 1),break;end
        [gcmma,xnew] = gcmma.innerUpdate(xnew,fxnew,gxnew, x,fx, dfdx, gx, dgdx, xmin,xmax);
        fxnew = xnew(1)^2 + xnew(2)^2 + xnew(3)^2;
        gxnew = [(xnew(1)-5)^2 + (xnew(2)-2)^2 + (xnew(3)-1)^2 - 9
            (xnew(1)-3)^2 + (xnew(2)-4)^2 + (xnew(3)-3)^2 - 9];

        [conserv] = gcmma.conCheck(fxnew,gxnew);
    end
    change = max(abs(x - xnew));
    x = xnew;
    disp(['iter: ' num2str(iter) ', change:' num2str(change) ', f(x): ', num2str(fx) ', x: ', num2str(x')])
end