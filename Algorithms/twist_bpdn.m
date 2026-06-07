function [v, cost] = twist_bpdn(H, u, lambda, maxIter, tol)
% TwIST for solving:
%   min_v (1/2)||Hv - u||_2^2 + lambda ||v||_1
%
% Inputs:
%   H        - forward model matrix
%   u        - measurements
%   lambda   - L1 regularisation parameter
%   maxIter  - maximum iterations
%   tol      - stopping tolerance
%
% Outputs:
%   v        - reconstructed solution
%   cost     - cost function history

[m, n] = size(H);

% Initialisation
v_prev = zeros(n,1);
v = v_prev;

% TwIST parameters (typical stable choice)
alpha = 1.0;
beta  = 0.5;

cost = zeros(maxIter,1);

for k = 1:maxIter
    
    % Gradient step
    grad = H'*(u - H*v);
    
    % Shrinkage step
    w = soft_thresh(v + grad, lambda);
    
    % TwIST update (two-step memory)
    v_new = (1 - alpha)*v_prev + (alpha - beta)*v + beta*w;
    
    % Cost function
    cost(k) = 0.5*norm(H*v_new - u)^2 + lambda*norm(v_new,1);
    
    % Convergence
    if k > 1 && abs(cost(k) - cost(k-1)) < tol
        cost = cost(1:k);
        break;
    end
    
    % Update variables
    v_prev = v;
    v = v_new;
end
end

% Soft threshold
function x = soft_thresh(z, thresh)
    x = sign(z).*max(abs(z) - thresh, 0);
end
