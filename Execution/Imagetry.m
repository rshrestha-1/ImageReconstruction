%% LSQR

maxit = 50;
tol = 1e-6;

v = lsqr(H,u,tol,maxit);

%% Sweep LSQR

for k = [5 10 15 20 30 40 50]
    v_temp = lsqr(H,u,tol,k);
    % Compute contrast and watch for noise and contrast degredation
end

%% L-curve for Optimal LSQR

% Inputs:
% H         : your system matrix (MxN)
% u         : measurement vector (Mx1)
% maxIter   : maximum iterations to test (e.g., 50)
% tol       : LSQR tolerance (e.g., 1e-6)

maxIter = 50;
tol = 1e-6;

residuals = zeros(maxIter,1);
solutionNorm = zeros(maxIter,1);

% Preallocate v
v = zeros(size(H,2),1);

for k = 1:maxIter
    % Run LSQR for k iterations
    [v_lsqr,~,~,~] = lsqr(H, u, tol, k);
    
    % Compute residual norm ||Hv - u||_2
    residuals(k) = norm(H*v_lsqr - u);
    
    % Compute solution norm ||v||_2
    solutionNorm(k) = norm(v_lsqr);
end

% Plot the L-curve
figure;
loglog(residuals, solutionNorm, '-o','LineWidth',2);
xlabel('Residual norm ||Hv - u||_2');
ylabel('Solution norm ||v||_2');
title('L-Curve for LSQR Reconstruction');
grid on;

% Annotate all iterations
for k = 1:maxIter
    text(residuals(k), solutionNorm(k), sprintf('%d',k), ...
        'FontSize',8, 'Color','b', 'HorizontalAlignment','left', 'VerticalAlignment','bottom');
end

% Highlight corner roughly
[~,cornerIdx] = max(diff(log10(residuals))./diff(log10(solutionNorm)));
hold on;
plot(residuals(cornerIdx), solutionNorm(cornerIdx), 'r*', 'MarkerSize',12);
text(residuals(cornerIdx), solutionNorm(cornerIdx), ...
    sprintf('  Optimal iteration ~ %d', cornerIdx), 'FontSize',12, 'Color','r');



