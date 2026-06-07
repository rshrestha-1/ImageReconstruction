function [v, cost] = twist_elastic(H, u, lambda1, lambda2, maxIter, tol)
% Elastic Net with TwIST:
%   min_v (1/2)||Hv - u||_2^2 + lambda1||v||_1 + lambda2||v||_2^2
%
% Inputs:
%   H        - forward model matrix
%   u        - measurements
%   lambda1   - L1 regularisation parameter
%   lambda2   - L2 regularisation parameter
%   maxIter  - maximum iterations
%   tol      - stopping tolerance
%
% Outputs:
%   v        - reconstructed solution
%   cost     - cost function history

[m, n] = size(H);

v_prev = zeros(n,1);
v = v_prev;

alpha = 1.0;
beta  = 0.5;

cost = zeros(maxIter,1);

for k = 1:maxIter
    
    % Gradient including L2 term
    grad = H'*(u - H*v) - 2*lambda2*v;
    
    % Shrinkage
    w = soft_thresh(v + grad, lambda1);
    
    % TwIST update
    v_new = (1 - alpha)*v_prev + (alpha - beta)*v + beta*w;
    
    % Cost
    cost(k) = 0.5*norm(H*v_new - u)^2 + lambda1*norm(v_new,1) + lambda2*norm(v_new,2)^2;
    
    if k > 1 && abs(cost(k) - cost(k-1)) < tol
        cost = cost(1:k);
        break;
    end
    
    v_prev = v;
    v = v_new;
end
end

function x = soft_thresh(z, thresh)
    x = sign(z).*max(abs(z) - thresh, 0);
end
