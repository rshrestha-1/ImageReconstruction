%% Carotid Phantom Simulation - assumes H and Hn in workspace
Nx = 19; Ny = 19; Nz = 11;

% Define range of voxel map
x_range = linspace(-2.2e-3, 2.2e-3, Nx); % m
y_range = linspace(-2.2e-3, 2.2e-3, Ny); % m
z_range = linspace(11.24e-3, 12.44e-3, Nz); % m

scenarios = ["healthy","lipid","fibrous","calcified","mixed"];

% Change accordingly with scenarios from 1 to 5
[v_true_3D2, labels] = carotid_phantom(Nx,Ny,Nz, scenarios(5));

%% Construct forwards model
v_true = v_true_3D2(:); % size = (Nx*Ny*Nz, 1)
u = Hn * v_true;

%% Simulate noise
SNR_dB = 20;

signal_power = mean(u.^2);
noise_power = signal_power / (10^(SNR_dB/10));

noise = sqrt(noise_power) * randn(size(u));

u = u + noise;
%% LSQR
tol = 1e-6;
maxit = 15;

v = lsqr(Hn,u,tol,maxit);
v_rec_3D2 = reshape(v,Nx,Ny,Nz);
%% FISTA
lambda_max = max(abs(H' * u));
lambda = 0.01 * lambda_max;

%lambdas = [0.1, 0.05, 0.01, 0.005] * lambda_max;

maxIter = 100;
tol = 1e-6;

[v, cost_fista] = fista_lasso(Hn, u, lambda, maxIter, tol);
v_rec_3D = reshape(v,Nx,Ny,Nz);
%% TwIST
% Normalise Hn
Hn_norm = Hn / norm(Hn);
u_n = u / norm(Hn);
lambda_max = max(abs(Hn_norm' * u));
lambda = 0.01 * lambda_max;

%lambdas = [0.1, 0.05, 0.01, 0.005] * lambda_max;

maxIter = 100;
tol = 1e-6;

[v, cost_twist] = twist_bpdn(Hn_norm, u_n, lambda, maxIter, tol);
v_rec_3D = reshape(v,Nx,Ny,Nz);
%% FISTA Elastic Net
lambda_max = max(abs(Hn' * u));
lambda1 = 0.1 * lambda_max;
lambda2 = 0.5 * lambda1;

maxIter = 300;
tol = 1e-6;

[v, cost_en_fista] = fista_elastic(Hn, u, lambda1, lambda2, maxIter, tol);
v_rec_3D = reshape(v,Nx,Ny,Nz);
%% TwIST Elastic Net
% Normalise Hn
Hn_norm = Hn / norm(Hn);
u_n = u / norm(Hn);
lambda_max = max(abs(Hn_norm' * u));
lambda1 = 0.025 * lambda_max;
lambda2 = 0.1 * lambda1;

maxIter = 400;
tol = 1e-6;

[v, cost_en_twist] = twist_elastic(Hn_norm, u_n, lambda1, lambda2, maxIter, tol);
v_rec_3D = reshape(v,Nx,Ny,Nz);

%% Reconstruct the images

% Envelope
V_env = abs(v_rec_3D);

% Normalise
V_env = V_env / max(V_env(:));

% Log compression (B-mode style)
dynamic_range_dB = 20;   % typical ultrasound range
V_dB = 20*log10(V_env + eps);  % avoid log(0)

% Clip dynamic range
V_dB(V_dB < -dynamic_range_dB) = -dynamic_range_dB;

% Plot all depth slices

figure;

for k = 1:Nz
    
    subplot(ceil(Nz/4), 4, k)  % 4 slices per row
    
    imagesc(x_range*1e3, y_range*1e3, squeeze(V_dB(:,:,k))')
    
    axis image
    set(gca,'YDir','normal')   % standard Cartesian
    colormap(gray)
    caxis([-dynamic_range_dB 0])
    
    title(sprintf('z = %.2f mm, at slice = %d', z_range(k)*1e3, k))
    
    if k == 1
        xlabel('x (mm)')
        ylabel('y (mm)')
    end
end

sgtitle('3D LSQR Reconstruction (Log Compressed B-mode)')
colorbar

% MIP (3 projections)

eps = 1e-6;
dynamic_range_dB = 10;
vol_mag = abs(v_true_3D);

% MIP along z (XY view)
% For each (x,y), take the maximum value along the z direction
mip_xy = max(vol_mag, [], 3);
mip_xy_db = 20*log10(mip_xy / max(mip_xy(:)) + eps);

% MIP along y (XZ view)
% For each (x,z), take the maximum value along the y direction
mip_xz = squeeze(max(vol_mag, [], 2));
mip_xz_db = 20*log10(mip_xz / max(mip_xz(:)) + eps);

% MIP along x (YZ view)
% For each (y,z), take the maximum value along the x direction
mip_yz = squeeze(max(vol_mag, [], 1));
mip_yz_db = 20*log10(mip_yz / max(mip_yz(:)) + eps);

% plot by voxel grid
figure

% XY projection (looking down the z-axis)
subplot(1,3,1)
imagesc(mip_xy_db.')
axis image
colormap gray
set(gca, 'YDir', 'normal');
caxis([-dynamic_range_dB 0])
colorbar
title('MIP - XY view (projection along Z)')
xlabel('X')
ylabel('Y')

% XZ projection (looking along the y-axis)
subplot(1,3,2)
imagesc(mip_xz_db)
axis image
colormap gray
caxis([-dynamic_range_dB 0])
colorbar
title('MIP - XZ view (projection along Y)')
xlabel('Z')
ylabel('X')

% YZ projection (looking along the x-axis)
subplot(1,3,3)
imagesc(mip_yz_db)
axis image
colormap gray
set(gca, 'YDir', 'normal');
caxis([-dynamic_range_dB 0])
colorbar
title('MIP - YZ view (projection along X)')
xlabel('Z')
ylabel('Y')

sgtitle('Maximum Intensity Projections of 3D Reconstruction by Voxel Grid')

% plot by dimension grid
figure

% XY projection (looking down the z-axis)
subplot(1,3,1)
imagesc(x_range*1e3, y_range*1e3, mip_xy_db.')
axis image
colormap gray
set(gca, 'YDir', 'normal');
caxis([-dynamic_range_dB 0])
colorbar
title('MIP - XY view (projection along Z)')
xlabel('X (mm)')
ylabel('Y (mm)')

% XZ projection (looking along the y-axis)
subplot(1,3,2)
imagesc(z_range*1e3, x_range*1e3, mip_xz_db)
axis image
colormap gray
caxis([-dynamic_range_dB 0])
colorbar
title('MIP - XZ view (projection along Y)')
xlabel('Z (mm)')
ylabel('X (mm)')

% YZ projection (looking along the x-axis)
subplot(1,3,3)
imagesc(z_range*1e3, y_range*1e3, mip_yz_db)
axis image
colormap gray
set(gca, 'YDir', 'normal');
caxis([-dynamic_range_dB 0])
colorbar
title('MIP - YZ view (projection along X)')
xlabel('Z (mm)')
ylabel('Y (mm)')

sgtitle('Maximum Intensity Projections of 3D Reconstruction by Dimension Grid')

% Compute metrics
results = compute_all_metrics(v_true_3D, v_rec_3D, labels);

fprintf('MSE = %d\nPSNR = %d\nSSIM = %d\nSNR = %d\nCNR = %d\nLipidContrast = %d\nCalcifiedContrast = %d\n', results.MSE, results.PSNR,results.SSIM,results.SNR,results.CNR,results.LipidContrast,results.CalcifiedContrast);