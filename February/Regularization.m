%% 
%% Example: L1 regularization
% Suppose H_f is m x n, u is m x 1

lambda_vals = logspace(-6, -1, 10);   % smaller lambdas if previous were too large
[B, FitInfo] = lasso(H_f, u, 'Lambda', lambda_vals);

%% Check sizes
disp(size(B))          % should be [n x numLambda]
disp(length(lambda_vals)) % should equal numLambda

%% L1 PLOT
figure;
semilogx(lambda_vals, B', 'LineWidth', 1.5);  % X = lambda, Y = v values
xlabel('Lambda (log scale)');
ylabel('v values');
title('Effect of L1 regularization on v');
grid on;
legend(arrayfun(@(x) sprintf('v%d',x), 1:size(H_f,2), 'UniformOutput', false), 'Location','best');
%% LSQR
tol = 1e-6;       
maxit = 130;      
[v, flag, relres, iter, resvec] = lsqr(H_f, u, tol, maxit);

figure; plot(resvec,'LineWidth',1.5);
xlabel('Iteration');
ylabel('Residual ||Hv - u||_2');
title('LSQR convergence');
grid on;
%% 
%% TWiST Visualization
H = H_f;
v = zeros(size(H,2),1);      
lambda = 0.1;                
alpha = 1 / (norm(H)^2);     
maxIter = 100;

residuals = zeros(maxIter,1);       
v_history = zeros(size(H,2), maxIter); 

for k = 1:maxIter
    % Gradient step (reduce L2 error)
    grad = H'*(H*v - u);
    v = v - alpha*grad;
    
    % Soft-thresholding (L1)
    v = sign(v).*max(abs(v)-lambda*alpha, 0);
    
    % Save data for plotting
    residuals(k) = norm(H*v - u);     
    v_history(:,k) = v;              
end

%% Plot residual per iteration
figure;
plot(1:maxIter, residuals);
xlabel('Iteration');
ylabel('Residual ||Hv - u||_2');
title('TWiST: Residual vs Iteration');
grid on;

%% Plot v components vs iteration
figure;
plot(1:maxIter, v_history', 'LineWidth',1.5);
xlabel('Iteration');
ylabel('v values');
title('TWiST: Evolution of v components');
grid on;
legend(arrayfun(@(x) sprintf('v%d',x), 1:size(H,2), 'UniformOutput', false));
