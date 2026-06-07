function [v, cost] = fista_lasso(H, u, lambda, maxIter, tol)
% FISTA for solving:
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

% Initialisation
[m, n] = size(H);
v = zeros(n,1);
y = v;
t = 1;

% Lipschitz constant L = ||H^T H||_2
L = norm(H)^2;

cost = zeros(maxIter,1);

for k = 1:maxIter
    
    % Gradient of data fidelity term
    grad = H'*(H*y - u);
    
    % Gradient descent + soft thresholding
    v_new = soft_thresh(y - (1/L)*grad, lambda/L);
    
    % Nesterov momentum update
    t_new = (1 + sqrt(1 + 4*t^2))/2;
    y = v_new + ((t - 1)/t_new)*(v_new - v);
    
    % Update variables
    v = v_new;
    t = t_new;
    
    % Cost function
    cost(k) = 0.5*norm(H*v - u)^2 + lambda*norm(v,1);
    
    % Convergence check
    if k > 1 && abs(cost(k) - cost(k-1)) < tol
        cost = cost(1:k);
        break;
    end
end
end

% Soft thresholding operator
function x = soft_thresh(z, thresh)
    x = sign(z).*max(abs(z) - thresh, 0);
end