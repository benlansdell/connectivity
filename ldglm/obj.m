function y = obj(x, S, f, Z, A, D, H, rho, X)
    
    foo = S.*log(f(x));
    foo(isnan(foo)) = 0; % 0*log(0) = 0 by convention
    y = sum(sum(f(x) - foo)) + sum(sum(Z.*A(x-D*H))) + rho/2*norm(A(x-D*H)-X,'fro')^2;
    
end