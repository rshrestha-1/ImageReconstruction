close all;
clear all;
clc;

load('H_5.mat');
load('udata_3points.mat');
u = udata_3points;
H_original = H_5;
rows = size(H_5,1);
%% Grid
Nx = 19;
Ny = 19;
Nz = 11;

N_voxels = Nx * Ny * Nz; % Total number of voxels

% Define range of voxel map
x_range = linspace(-2.2e-3, 2.2e-3, Nx); % mm
y_range = linspace(-2.2e-3, 2.2e-3, Ny); % mm 
z_range = linspace(11.24e-3, 12.44e-3, Nz); % mm

[x_grid, y_grid, z_grid] = ndgrid(x_range, y_range, z_range);

voxel_grid = reshape(H_original, [rows, Nx, Ny, Nz]);

%% Parameters
D_values = [1 2 3 4 5];
max_iter_values = [30 30 30 30];

Ne = 5;
Nt = rows / Ne;
tol = 1e-8;
 

%%  STORAGE
results_mse  = [];
results_mae  = [];
results_ssim = [];
row = 1;

%%  MAIN LOOP
for D = D_values
    
    % Downsample
    x_ds = x_range(1:D:end);
    y_ds = y_range(1:D:end);
    z_ds = z_range(:);

    voxel_ds = voxel_grid(:, 1:D:end, 1:D:end, :);

    H_5_down = zeros(rows, Nx*Ny*Nz);

    for i = 1:rows
        Vtemp = squeeze(voxel_ds(i,:,:,:));
        Vq = interpn(x_ds, y_ds, z_ds, Vtemp, x_grid, y_grid, z_grid, 'linear');
        H_5_down(i,:) = Vq(:)';
    end

    % Convolution
    Hf_3D = reshape(H_5_down, [Nt, Ne, N_voxels]);
    H_conv = zeros(Ne*Nt, N_voxels);
    

    for v_ind = 1:N_voxels
        for e = 1:Ne
            h_tx = Hf_3D(:, e, v_ind);
            h_echo = conv(h_tx, h_tx);
            h_echo = h_echo(1:Nt);
            h_echo = h_echo / max(abs(h_echo) + eps);

            row_start = (e-1)*Nt + 1;
            row_end   = e*Nt;

            H_conv(row_start:row_end, v_ind) = h_echo;
        end
    end

    for max_iter = max_iter_values
        if D == 1
            v_ref = lsqr(H_conv, u, tol, max_iter);
            V_ref = reshape(abs(v_ref), Nx, Ny, Nz);

        end
        
        % Reconstruction
        v_est = lsqr(H_conv, u, tol, max_iter);
        V_est = reshape(abs(v_est), Nx, Ny, Nz);
        


        % MSE (mean + std)
        mse_vals = (abs(v_est) - abs(v_ref)).^2;
        mse_val = mean(mse_vals);
        mse_std = std(mse_vals);
        
        % MAE (mean + std)
        mae_vals = abs(abs(v_est) - abs(v_ref));
        mae_val = mean(mae_vals);
        mae_std = std(mae_vals);
        
        % SSIM (mean + std across slices)
        V_est_n = V_est / (max(V_est(:)) + eps);
        V_ref_n = V_ref / (max(V_ref(:)) + eps);
        
        ssim_vals = zeros(1, Nz);
        for k = 1:Nz
            ssim_vals(k) = ssim(V_est_n(:,:,k), V_ref_n(:,:,k));
        end
        
        ssim_mean = mean(ssim_vals);
        ssim_std  = std(ssim_vals);

        % Store
        results_mse(row,:)  = [D, max_iter, mse_val, mse_std];
        results_mae(row,:)  = [D, max_iter, mae_val, mae_std];
        results_ssim(row,:) = [D, max_iter, ssim_mean, ssim_std];

        row = row + 1;
    end
end

%% TABLES

%% ===== CONVERT TO MATRIX FORMAT =====

D_vals = unique(results_mse(:,1));
iter_vals = unique(results_mse(:,2));

nD = length(D_vals);
nI = length(iter_vals);

MSE_mean  = zeros(nD, nI);
MSE_std   = zeros(nD, nI);

MAE_mean  = zeros(nD, nI);
MAE_std   = zeros(nD, nI);

SSIM_mean = zeros(nD, nI);
SSIM_std  = zeros(nD, nI);

for i = 1:size(results_mse,1)

    D    = results_mse(i,1);
    iter = results_mse(i,2);

    row = find(D_vals == D);
    col = find(iter_vals == iter);

    MSE_mean(row,col)  = results_mse(i,3);
    MSE_std(row,col)   = results_mse(i,4);
    
    MAE_mean(row,col)  = results_mae(i,3);
    MAE_std(row,col)   = results_mae(i,4);
    
    SSIM_mean(row,col) = results_ssim(i,3);
    SSIM_std(row,col)  = results_ssim(i,4);

end
%% 

formatCell = @(mean_mat, std_mat) ...
    arrayfun(@(m,s) sprintf('%.3e ± %.3e', m, s), ...
    mean_mat, std_mat, 'UniformOutput', false);

MSE_disp  = formatCell(MSE_mean, MSE_std);
MAE_disp  = formatCell(MAE_mean, MAE_std);
SSIM_disp = formatCell(SSIM_mean, SSIM_std);


table_mse = cell2table(MSE_disp, ...
    'VariableNames', compose('Iter_%d', iter_vals), ...
    'RowNames', compose('D_%d', D_vals));

table_mae = cell2table(MAE_disp, ...
    'VariableNames', compose('Iter_%d', iter_vals), ...
    'RowNames', compose('D_%d', D_vals));

table_ssim = cell2table(SSIM_disp, ...
    'VariableNames', compose('Iter_%d', iter_vals), ...
    'RowNames', compose('D_%d', D_vals));
%% 
disp('--- MSE (rows = D, cols = iterations) ---')
disp(table_mse)

disp('--- MAE (rows = D, cols = iterations) ---')
disp(table_mae)

disp('--- SSIM (rows = D, cols = iterations) ---')
disp(table_ssim)
