a = 2;
b = 0.5;

xmin = -3;
xmax = 2;
xmid = 0.5*(xmin+xmax);

X = xmin:0.01:xmax;
plot(X,exp(a*(X-b))-a*(X-b) - 1,'k','LineWidth',3);
xlabel('x'); ylabel('f(x)')
set(gca,'FontSize',18)
change = 1;
for iter = 1:50
    if change < 1e-5,break;end
    xold = xmid;
    
    fx = exp(a*(xmid-b))-a*(xmid-b) - 1;
    dfdx = a*exp(a*(xmid-b))-a;
    
    if dfdx < 0,xmin = xmid;
    elseif dfdx > 0,xmax = xmid;
    end 
    xmid = 0.5*(xmin+xmax);
    change = max(abs(xmid - xold));     
    disp(['iter: ' num2str(iter) ', change:' num2str(change) ', x: ', num2str(xmid') ', f(x): ', num2str(fx)])
end

