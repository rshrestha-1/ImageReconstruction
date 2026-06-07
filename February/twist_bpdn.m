function v = twist_bpdn(A, AT, u, lambda, maxit,tau)


% TwIST parameters
alpha = 1;
beta  = 0.5;

% Initialization
v_prev = zeros(size(AT(u)));
v      = v_prev;

soft = @(x,t) sign(x).*max(abs(x)-t,0);

for k = 1:maxit
    
    grad = AT(A(v) - u);
    v_new = soft(v - tau*grad, tau*lambda);
    
    if k > 1
        v_new = alpha*v_new + beta*(v_new - v_prev);
    end
    
    v_prev = v;
    v = v_new;
end

end
