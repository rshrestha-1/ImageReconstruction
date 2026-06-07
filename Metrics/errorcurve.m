%% eror curve
lambda1_vals = linspace(0.001, 0.1, 20) * max(abs(v_l2));
errors = zeros(size(lambda1_vals));

for i = 1:length(lambda1_vals)
    
    lambda1 = lambda1_vals(i);
    
    % Apply L1 shrinkage
    v_test = sign(v_l2) .* max(abs(v_l2) - lambda1, 0);
    
    % Find reconstructed voxel
    [~, idx] = max(abs(v_test));
    
    % Compute error
    errors(i) = abs(idx - 1625);
end

plot(lambda1_vals, errors, 'LineWidth', 2);
xlabel('\lambda_1');
ylabel('Localization Error');
title('Error vs L1 Regularisation');
grid on;

[~, best_idx] = min(errors);
best_lambda1 = lambda1_vals(best_idx);
%% SNR plot
peak = max(abs(v_test));
noise = mean(abs(v_test));
SNR = peak / noise;

%Plot SNR
%% heatmap lambda 1 and 2
lambda1_vals = linspace(0.01, 0.1, 15);   % L1
lambda2_vals = logspace(-5, -2, 15);      % L2 (log scale is better)

error_map = zeros(length(lambda2_vals), length(lambda1_vals));

for i = 1:length(lambda2_vals)
    for j = 1:length(lambda1_vals)
        
        lambda2 = lambda2_vals(i);
        lambda1 = lambda1_vals(j);
        
        % --- L2 (LSQR with Tikhonov) ---
        A = [H; sqrt(lambda2)*eye(size(H,2))];
        b_aug = [b; zeros(size(H,2),1)];
        
        v_l2 = lsqr(A, b_aug, 1e-6, 30);
        
        % --- L1 (soft threshold) ---
        v_sparse = sign(v_l2) .* max(abs(v_l2) - lambda1*max(abs(v_l2)), 0);
        
        % --- Find peak voxel ---
        [~, idx] = max(abs(v_sparse));
        
        % --- Compute error ---
        error_map(i,j) = abs(idx - 1625);
    end
end

imagesc(lambda1_vals, log10(lambda2_vals), error_map);

colorbar;
xlabel('\lambda_1 (L1 sparsity)');
ylabel('log_{10}(\lambda_2) (L2 regularisation)');
title('Reconstruction Error Heatmap');

set(gca, 'YDir', 'normal');
[min_err, idx] = min(error_map(:));
[i_opt, j_opt] = ind2sub(size(error_map), idx);

best_lambda2 = lambda2_vals(i_opt);
best_lambda1 = lambda1_vals(j_opt);
%% Normalise lambda
lambda1_vals = linspace(0.01, 0.1, 15);
lambda1_scaled = lambda1 * max(abs(v_l2));
% contour lines
hold on;
contour(lambda1_vals, log10(lambda2_vals), error_map, 10, 'k');
%% CRLB fischer
sigma2 = var(noise);  % estimate noise variance

FIM = (1/sigma2) * (H' * H);
CRLB = sigma2 * inv(H' * H);