%% Clear and load
close all;
clear all;
clc;

load('H_5.mat');           % Original voxel data
load('udata_3points.mat_u.mat'); % Other data (not used here yet)

%% Original grid size
Nx = 19;
Ny = 19;
Nz = 11;
N_voxels = Nx * Ny * Nz; 

% voxel coordinates
x_range = linspace(-2.2e-3, 2.2e-3, Nx);
y_range = linspace(-2.2e-3, 2.2e-3, Ny);
z_range = linspace(11.24e-3, 12.44e-3, Nz);

% Create 3D grid
[x_grid, y_grid, z_grid] = ndgrid(x_range, y_range, z_range);

voxel_grid = reshape(H_5, [4575, Nx, Ny, Nz]);

%% Downsample parameters
D = 3; % downsampling factor

% Downsample voxel coordinates
x_ds = x_range(1:D:end);
y_ds = y_range(1:D:end);
z_ds = z_range(1:D:end);

% Downsample the voxel data
voxel_ds = voxel_grid(:, 1:D:end, 1:D:end, 1:D:end);
[new_Nx, new_Ny, new_Nz] = size(voxel_ds(1,:,:,:));

new_Nx = squeeze(new_Nx);
new_Ny = squeeze(new_Ny);
new_Nz = squeeze(new_Nz);

%% Interpolate each volume back to original grid using vector mode
H_5_down = zeros(4575, Nx*Ny*Nz);

[xq_grid, yq_grid, zq_grid] = ndgrid(x_range, y_range, z_range);

for i = 1:4575
    V = squeeze(voxel_ds(i,:,:,:)); % downsampled volume
   
    Vq = interp3(x_ds, y_ds, z_ds, V, xq_grid, yq_grid, zq_grid, 'linear'); 
    
    H_5_down(i,:) = Vq(:)';
end
%% 

Ne = 5;
Nt = size(H_5_down,1)/Ne;


Hf_3D = reshape(H_5_down, [Nt, Ne, N_voxels]);

H_conv = zeros(Ne*Nt, N_voxels);

for v = 1:N_voxels
    
    for e = 1:Ne
      
        h_tx = Hf_3D(:, e, v);
        
        % autoconvolution
        h_echo_full = conv(h_tx, h_tx);
        
        % truncate
        h_echo = h_echo_full(1:Nt);
        
        % Normalisation
        h_echo = h_echo / max(abs(h_echo) + eps);
        
        % put into initialised H
        row_start = (e-1)*Nt + 1;
        row_end   = e*Nt;
        
        H_conv(row_start:row_end, v) = h_echo;
        
    end
end

%% 

max_iter = 15;
tol = 1e-8;
v = lsqr(H_conv,udata_3points.mat, tol, max_iter);
%% 
V = reshape(abs(v), Nx, Ny, Nz);

V_env_1 = abs(V);

% --- Normalise
V_env_1 = V_env_1 / max(V_env_1(:));

% --- Log compression (B-mode style)
dynamic_range_dB = 15;  
V_dB = 20*log10(V_env_1 + eps);  

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
    datacursormode on
    title(sprintf('z = %.2f mm, at slice = %d', z_range(k)*1e3, k))
    
    if k == 1

        xlabel('x (mm)')
        ylabel('y (mm)')
    end
end
sgtitle(['Reconstruction with ' num2str(max_iter) ' iterations and ' num2str(dynamic_range_dB) 'dB' ])
%%