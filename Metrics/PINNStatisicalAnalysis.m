%% PINNs Statistical Analysis
H_pred = H_time;
H_true = H_exp_new;
%% Plot training loss convergence

figure;

% Plot loss history
plot(loss_total, 'LineWidth', 1.5);

% Use logarithmic scale for y-axis
set(gca, 'YScale', 'log');

% Labels and title
xlabel('Epoch');
ylabel('Loss');
title('Training Convergence');

% Grid for clarity
grid on;
set(gca, 'FontSize', 12);

%% Plot multiple loss convergence

figure; hold on;

% Plot each loss
plot(loss_total, 'LineWidth', 1.8);
plot(loss_data, '--', 'LineWidth', 1.5);
plot(loss_physics, ':', 'LineWidth', 1.5);

% Log scale for better visibility
set(gca, 'YScale', 'log');

% Labels
xlabel('Epoch');
ylabel('Loss');

% Title
title('PINN Training Convergence');

% Legend
legend('Total Loss', 'Data Loss', 'Physics Loss', 'Location', 'best');

% Grid
grid on;

% Formatting
set(gca, 'FontSize', 12);
axis tight;

%% Compute metrics: MSE, RMSE, NRMSE, Correlation Coefficient

% Flatten matrices into vectors
true = H_true(:);
pred = H_pred(:);

% Mean Squared Error (MSE)
mse = mean(abs(true - pred).^2);

% Root Mean Squared Error (RMSE)
rmse = sqrt(mse);

% Normalised RMSE
nrmse = rmse / max(abs(true));

% Correlation Coefficient

% corrcoef returns a 2x2 matrix, take off-diagonal
C = corrcoef(abs(true), abs(pred));
corr_val = C(1,2);

% Display results
fprintf('MSE: %e\n', mse);
fprintf('RMSE: %e\n', rmse);
fprintf('NRMSE: %e\n', nrmse);
fprintf('Correlation: %.4f\n', corr_val);

%% Set freqs for interpolation

freqs = [322265.625, 377197.265625, 373535.15625, 380859.375, 25634.765625, 318603.515625, ...
    292968.75, 289306.640625, 296630.859375, 300292.96875, 314941.40625, 303955.078125, ...
    311279.296875, 307617.1875, 21972.65625, 18310.546875, 14648.4375, 10986.328125, 7324.21875, 3662.109375];

% Reshape H
num_sensors = 58;
Nt = 8192;
Nx = 20;
Ny = 20;
Nz = 20;

dt = 1/30e6;

H_time_reshaped = reshape(H_time, [num_sensors, Nt, Nx*Ny*Nz]);

% Convert to frequency domain
H_fft = fft(H_time_reshaped, [], 2);
H_fft = reshape(H_fft, [num_sensors, Nt, Nx, Ny, Nz]);

% Build full frequency axis
fs = 1/dt; % you must know dt
freqs_full = (0:Nt-1)*(fs/Nt);
freqs_full(freqs_full > fs/2) = freqs_full(freqs_full > fs/2) - fs;

% Find indices of your selected frequencies
idx = zeros(length(freqs),1);

for i = 1:length(freqs)
    [~, idx(i)] = min(abs(freqs_full - freqs(i)));
end

% Sort frequencies and indices together
[freqs, sort_idx] = sort(freqs);
idx = idx(sort_idx);

% Extract those frequencies
H_selected = H_fft(:, idx, :, :, :);
%% Physics residual

c = 1500;  % wave speed

% Set voxel spacing
    dx = 0.0003; % 0.3mm

[num_sensors, K, Nx, Ny, Nz] = size(H_selected);

% Preallocate residual storage
residual_normalised_all = zeros(K,1);

% Loop frequencies

for i = 1:K

    % Get frequency
    f = freqs(i);

    % Compute wavenumber
    k = 2*pi*f / c;

    % Loop and average sensors
    res_temp = zeros(num_sensors,1);

    for s = 1:num_sensors
    
        u = squeeze(H_selected(s, i, :, :, :));
    
        % Laplacian 3D
        lap = ( ...
            circshift(u,[1,0,0]) + circshift(u,[-1,0,0]) + ...
            circshift(u,[0,1,0]) + circshift(u,[0,-1,0]) + ...
            circshift(u,[0,0,1]) + circshift(u,[0,0,-1]) ...
            - 6*u ) / dx^2;
    
        % Helmoltz Residual
        residual = lap + (k^2) * u;
    
        % Normalised Residual
        res_temp(s) = norm(residual(:)) / norm(u(:));
    
    end

    % Average across sensors
    residual_normalised_all(i) = mean(res_temp);

end

% Print Results

fprintf('Physics Residual Summary \n');
fprintf('Mean residual: %e\n', mean(residual_normalised_all));
fprintf('Min residual:  %e\n', min(residual_normalised_all));
fprintf('Max residual:  %e\n', max(residual_normalised_all));

% Plot Residual vs Frequency

figure;

plot(freqs, residual_normalised_all, 'LineWidth', 1.5);

xlabel('Frequency (Hz)');
ylabel('Normalised Residual');
title('Helmholtz Residual vs Frequency');

grid on;
set(gca, 'FontSize', 12);

