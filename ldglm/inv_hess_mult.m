function y = inv_hess_mult(H,x, T, rho, center, N)
    % For the Newton step of ADMM, we must multiply a vector by the inverse of
    % H + rho*A'*A, where H is the Hessian of the smooth penalty term (which in
    % our case is diagonal and therefore simple) and A is the linear operator
    % in the objective ||A(x)||_* + f(x). In our problem the linear operator is
    % mean-centering the matrix x, which is a symmetric and idempotent operator
    % that can be written out in matrix form for m-by-n matrices as:
    %
    % A = eye(m*n) - 1/n*kron(ones(n,1),eye(m))*kron(ones(n,1),eye(m))'
    %
    % and so (H + rho*A*A')^-1*x can be efficiently computed by taking
    % advantage of Woodbury's lemma.
    
    %T
    %rho 
    %center 

    if center % update this to include spike history term
        Hi = 1./(H + rho);
        foo = (sum(Hi.*x,2)./(ones(N,1)-sum(Hi,2)*rho/T));
        y = Hi.*x + rho/T*Hi.*(foo*ones(1,T));
    else
        y = x./(H + rho);
    end
end