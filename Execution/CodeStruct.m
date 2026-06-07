%% RESET
close all; clear all; clc
%% 1a. Set H
H = H_5;
%% 1b. Set u
u = udata_3points;
%% 2. Voxel Map
% Define number of voxels
Nx = 19; % number of lateral voxels
Ny = 19; % number of height/elevational voxels
Nz = 11; % number of depth voxels

N_voxels = Nx * Ny * Nz; % Total number of voxels

% Define range of voxel map
x_range = linspace(-2.2e-3, 2.2e-3, Nx); % m
y_range = linspace(-2.2e-3, 2.2e-3, Ny); % m
z_range = linspace(11.24e-3, 12.44e-3, Nz); % m

[x_grid, y_grid, z_grid] = ndgrid(x_range, y_range, z_range);

voxel_positions = [x_grid(:), y_grid(:), z_grid(:)];
%% 3a. Normalise H
H = H / max(abs(H(:)));
%% 3b. Normalise u
u = u / max(abs(u));
%% 4. Autoconvolution

Ne = 5; % number of elements
Nt = 4575/Ne; % number of time samples per element
Nv = Nx*Ny*Nz; % number of voxels

eps = 1e-6;

Hf_3D = reshape(H, [Nt, Ne, Nv]);
H = zeros(Ne*Nt, Nv);

for v = 1:Nv
    
    for e = 1:Ne
      
        h_tx = Hf_3D(:, e, v);
        
        % Autoconvolution
        h_echo_full = conv(h_tx, h_tx);
        
        % Truncate
        h_echo = h_echo_full(1:Nt);
        
        % Normalisation
        h_echo = h_echo / max(abs(h_echo) + eps);
        
        % Into initialised H
        row_start = (e-1)*Nt + 1;
        row_end   = e*Nt;
        
        H(row_start:row_end, v) = h_echo;
        
    end
    
end
%% Bandpass and u Gaussian filtering
fs = 30e6;

f_low = 10.39e6;
f_high = 13e6;

[b,a] = butter(4, [f_low f_high]/(fs/2), 'bandpass'); % normalised to Nyquist
u_filt = filtfilt(b,a,u);
u_filt = smoothdata(u_filt, 'gaussian', 5);
%u_filt = u_filt / max(abs(u_filt));
%% LSQR
tol = 1e-6;
maxit = 15;

v = lsqr(Hn,u_filt,tol,maxit);
volume = reshape(v,Nx,Ny,Nz);
%% FISTA
lambda_max = max(abs(Hn' * u));
lambda = 0.01 * lambda_max;

%lambdas = [0.1, 0.05, 0.01, 0.005] * lambda_max;

maxIter = 100;
tol = 1e-6;

[v, cost] = fista_lasso(Hn, u, lambda, maxIter, tol);
volume = reshape(v,Nx,Ny,Nz);
%% TwIST
% Normalise Hn
Hn_norm = Hn / norm(Hn);
u_n = u / norm(Hn);
lambda_max = max(abs(Hn_norm' * u));
lambda = 0.01 * lambda_max;

%lambdas = [0.1, 0.05, 0.01, 0.005] * lambda_max;

maxIter = 100;
tol = 1e-6;

[v, cost] = twist_bpdn(Hn_norm, u_n, lambda, maxIter, tol);
volume = reshape(v,Nx,Ny,Nz);
%% FISTA Elastic Net
lambda_max = max(abs(Hn' * u));
lambda1 = 0.1 * lambda_max;
lambda2 = 0.5 * slambda1;

maxIter = 300;
tol = 1e-6;

[v, cost] = fista_elastic(Hn, u, lambda1, lambda2, maxIter, tol);
volume = reshape(v,Nx,Ny,Nz);
%% TwIST Elastic Net
% Normalise Hn
Hn_norm = Hn / norm(Hn);
u_n = u / norm(Hn);
lambda_max = max(abs(Hn_norm' * u));
lambda1 = 0.01 * lambda_max;
lambda2 = 0.1 * lambda1;

maxIter = 200;
tol = 1e-6;

[v, cost] = twist_elastic(Hn_norm, u_n, lambda1, lambda2, maxIter, tol);
volume = reshape(v,Nx,Ny,Nz);

%% Avoid for Integrity to Recon

% Threshold
v(abs(v) < 0.05*max(abs(v))) = 0;

volume = reshape(v,Nx,Ny,Nz);

v_smooth = imgaussfilt3(volume, 1); % sigma = 1 voxel
v_med = medfilt3(volume, [3 3 3]);

%% B-mode (all XY slices)

% Envelope
V_env = abs(volume);

% Normalise
V_env = V_env / max(V_env(:));

% Log compression (B-mode style)
dynamic_range_dB = 15;   % typical ultrasound range
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

%% MIP (3 projections)

eps = 1e-6;
dynamic_range_dB = 10;
vol_mag = abs(volume);

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
%% Plot v
figure;
plot(v, 'LineWidth', 1.5);
xlabel('Index');
ylabel('Intensity');
title('Intensity Values');
grid on;
