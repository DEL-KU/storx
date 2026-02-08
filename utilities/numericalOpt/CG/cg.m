function u = cg(K,f,u,tol)
r = f - K*u;
fnorm = norm(f);
p = r;
k = 0;
res = norm(r)/fnorm;
while res > tol
    r_old = r;
    alpha = (r_old'*r_old)/(p'*K*p);
    r = r - alpha*K*p;
    u = u + alpha*p;
    beta = (r'*r)/(r_old'*r_old);
    p = r + beta*p;
    k = k+1;
    res = norm(r)/fnorm;
    disp(['iter: ' num2str(k) ', residual:' num2str(res) ', u: ', num2str(u')])
end

