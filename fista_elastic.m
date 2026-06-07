function [v, cost] = fista_elastic(H, u, lambda1, lambda2, maxIter, tol)
% Elastic Net via FISTA:
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

v = zeros(n,1);
y = v;
t = 1;

% Modified Lipschitz constant (includes L2 term)
L = norm(H)^2 + 2*lambda2;

cost = zeros(maxIter,1);

for k = 1:maxIter
    
    % Gradient including L2 term
    grad = H'*(H*y - u) + 2*lambda2*y;
    
    % Proximal step (L1 only)
    v_new = soft_thresh(y - (1/L)*grad, lambda1/L);
    
    % Nesterov momentum update
    t_new = (1 + sqrt(1 + 4*t^2))/2;
    y = v_new + ((t - 1)/t_new)*(v_new - v);
    
    v = v_new;
    t = t_new;
    
    % Cost
    cost(k) = 0.5*norm(H*v - u)^2 + lambda1*norm(v,1) + lambda2*norm(v,2)^2;
    
    if k > 1 && abs(cost(k) - cost(k-1)) < tol
        cost = cost(1:k);
        break;
    end
end
end

function x = soft_thresh(z, thresh)
    x = sign(z).*max(abs(z) - thresh, 0);
end