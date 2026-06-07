%% Function Handles
A = @(x)H*x;
AT = @(x)H'*x;

%% Verify Adjoint Consistency

x = randn(size(H,2),1);
y = randn(size(H,1),1);
left = A(x)'*y;
right = x'*AT(y);
disp(left-right)

%% LSQR

maxit = 15;
tol = 1e-6;

v = lsqr(A,u,tol,maxit);
