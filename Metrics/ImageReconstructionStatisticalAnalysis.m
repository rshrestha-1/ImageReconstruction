%% Reconstruction Algorithm Statistical Analysis

%% Cost function convergence

figure;
hold on; grid on; grid minor;

% Plot cost histories (log scale for convergence clarity)
h1 = plot(cost_en_fista, 'LineWidth', 2);
h2 = plot(cost_en_twist, 'LineWidth', 2);
h3 = plot(cost_fista, 'LineWidth', 2);
h4 = plot(cost_twist, 'LineWidth', 2);

% Log scale
set(gca, 'YScale', 'log');

% Labels
xlabel('Iterations', 'FontSize', 13);
ylabel('Cost Function Value', 'FontSize', 13);
title('Algorithm Convergence Comparison', 'FontSize', 14);

% Legend
legend([h1 h2 h3 h4], ...
    {'FISTA', 'TwIST', 'Elastic Net FISTA', 'Elastic Net TwIST'}, ...
    'Location', 'northeast', 'FontSize', 10);

% Aesthetic improvements
set(gca, 'FontSize', 12);
box on;
xlim([1 100]);
ylim([3e-3 1e0]);

% Markers
set(h1, 'Marker', 'o', 'MarkerIndices', 1:10:length(cost_fista));
set(h2, 'Marker', 's', 'MarkerIndices', 1:10:length(cost_twist));
set(h3, 'Marker', 'd', 'MarkerIndices', 1:10:length(cost_en_fista));
set(h4, 'Marker', '^', 'MarkerIndices', 1:10:length(cost_en_twist));

%% MSE

mse_val = mean((x(:) - v_rec(:)).^2)

%% SNR

signal_power = norm(x(:))^2;
noise_power  = norm(x(:) - v_rec(:))^2;

snr_val = 10 * log10(signal_power / noise_power)

%% NRMSE

nrmse = norm(x(:) - v_rec(:)) / norm(x(:))

%% Difference Image
scale_factor = 5; % scaling for visibility

% Compute absolute difference
diff_img = abs(v_true - v_rec) * scale_factor;

% Display side-by-side
figure;

subplot(1,3,1);
imagesc(v_true);
axis image; colormap gray; colorbar;
title('Ground Truth');

subplot(1,3,2);
imagesc(v_rec);
axis image; colormap gray; colorbar;
title('Reconstruction');

subplot(1,3,3);
imagesc(diff_img);
axis image; colormap gray; colorbar;
title(['Absolute Difference x' num2str(scale_factor)]);
%% 
scale_factor = 5;
slice_num = 5; % pick a slice
diff_img = abs(vol(:,:,slice_num) - volume(:,:,slice_num)) * scale_factor;

figure;
imagesc(diff_img);
axis image; colormap hot; colorbar;
title(['Absolute Difference (Slice ' num2str(slice_num) ') x' num2str(scale_factor)]);


