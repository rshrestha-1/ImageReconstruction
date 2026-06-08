close all; 
clear all; 
clc

load('H_exp_new.mat');
load('u_exp_new.mat');
%% 

H_original = H_exp_new;
u = u_exp_new;
%% ORIGINAL VOXEL GRID (USED BY H)

Nx = 6; 
Ny = 6;
Nz = 6; 

N_voxels = Nx * Ny * Nz; 

x_range = linspace(-2.5e-3, 2.5e-3, Nx); 
y_range = linspace(-2.5e-3, 2.5e-3, Ny); 
z_range = linspace(14e-3, 19e-3, Nz); 

[x_grid, y_grid, z_grid] = ndgrid(x_range, y_range, z_range);

%% PARAMETERS

Ne = 58;
Nt = size(H_original,1)/Ne;

%% RESHAPE H 

H_5D = reshape(H_original, [Nt, Ne, Nx, Ny, Nz]);

%% TAKE FRONT SLICE 

H_front = H_5D(:,:,:,:,1);   

%% ANGULAR SPECTRUM PARAMETERS

dx = 0.001;       
dy = 0.001;
      
lambda = 0.2e-3; 

%frequency index
k = 2*pi/lambda; 

% Fourier coordinates
fx = (-Nx/2:Nx/2-1)/(Nx*dx);
fy = (-Ny/2:Ny/2-1)/(Ny*dy);
[FX,FY] = meshgrid(fx, fy);

%Kz equation
kz = sqrt(k^2 - (2*pi*FX).^2 - (2*pi*FY).^2);

dz = 1e-3;  

%number of slices we want to propagate to
Nz_asp = 2;

z_idx = 5; % choose your reference plane
Hf_3D = reshape(H_original, [Nt, Ne, Nx, Ny, Nz]);

H_xy = Hf_3D(:,:,:,:,z_idx);  
% size: Nt × Ne × Nx × Ny

%% PROPAGATE H USING ANGULAR SPECTRUM

H_prop = zeros(Nt, Ne, Nx, Ny, Nz_asp);

for t = 1:Nt
    for e = 1:Ne
        
        slice0 = squeeze(H_xy(t,e,:,:)); % Nx × Ny
        
        F_slice = fftshift(fft2(slice0));
        
        for zi = 1:Nz_asp
            F_prop = F_slice .* exp(1i * kz * dz * zi);
            H_prop(t,e,:,:,zi) = ifft2(ifftshift(F_prop));
        end
        
    end
end

%% RESHAPE BACK

N_voxels_new = Nx*Ny*Nz_asp;

%% AUTOCONVOLUTION OF H 

H_conv = zeros(Ne*Nt, N_voxels_new);
Hf_3D = reshape(H_prop, [Nt, Ne, N_voxels_new]);


for v = 1:N_voxels_new
    
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

z_range_asp = z_range(5) + (1:Nz_asp)*dz;

%% LSQR RECONSTRUCTION
% col_norms = sqrt(sum(H_conv.^2,1));
% Hn = H_conv ./ col_norms;

%% 

maxIter = 15;
tol = 1e-8;

v = lsqr(H_conv, u, tol, maxIter);

%% RESHAPE RESULT

V = reshape(abs(v), Nx, Ny, Nz_asp);

%% DISPLAY RESULT

V_env_1 = abs(V);

% --- Normalise
V_env_1 = V_env_1 / max(V_env_1(:));

% --- Log compression (B-mode style)
dynamic_range_dB = 10;   % typical ultrasound range
V_dB = 20*log10(V_env_1 + eps);  % avoid log(0)

% Clip dynamic range
V_dB(V_dB < -dynamic_range_dB) = -dynamic_range_dB;

% Plot all depth slices

figure;

for k = 1:Nz_asp
    
    subplot(ceil(Nz_asp/4), 4, k)  % 4 slices per row
    
    imagesc(x_range*1e3, y_range*1e3, squeeze(V_dB(:,:,k))')
    
    axis image
    set(gca,'YDir','normal')   % standard Cartesian
    colormap(gray)
    caxis([-dynamic_range_dB 0])
    datacursormode on
    title(sprintf('z = %.2f mm, at slice = %d', z_range_asp(k)*1e3, k))
    
    if k == 1
        xlabel('x (mm)')
        ylabel('y (mm)')
    end
end
sgtitle(['Reconstruction with ' num2str(maxIter) ' iterations and ' num2str(dynamic_range_dB) 'dB (trial 3)' ])
%% Original 3D voxel volume

Nx_fine = 20;
Ny_fine = 20;
Nz_fine = 20;

x_fine = linspace(min(x_range), max(x_range), Nx_fine);
y_fine = linspace(min(y_range), max(y_range), Ny_fine);
z_fine = linspace(min(z_range_asp), max(z_range_asp), Nz_fine);
[x_grid, y_grid, z_grid] = ndgrid(x_range, y_range, z_range_asp);
[xq, yq, zq] = ndgrid(x_fine, y_fine, z_fine);

V_angular = interpn(x_grid, y_grid, z_grid, V, xq, yq, zq, 'linear', 0);
%% 

V_env_1 = abs(V_angular);

% --- Normalise
V_env_1 = V_env_1 / max(V_env_1(:));

% --- Log compression (B-mode style)
dynamic_range_dB = 15;   % typical ultrasound range
V_dB = 20*log10(V_env_1 + eps);  % avoid log(0)

% Clip dynamic range
V_dB(V_dB < -dynamic_range_dB) = -dynamic_range_dB;

% Plot all depth slices

figure;

for k = 1:Nz_asp
    
    subplot(ceil(Nz_asp/4), 4, k)  % 4 slices per row
    
    imagesc(x_fine*1e3, y_fine*1e3, squeeze(V_dB(:,:,k))')
    
    axis image
    set(gca,'YDir','normal')   % standard Cartesian
    colormap(gray)
    caxis([-dynamic_range_dB 0])
    datacursormode on
    title(sprintf('z = %.2f mm, at slice = %d', z_range_asp(k)*1e3, k))
    
    if k == 1
        xlabel('x (mm)')
        ylabel('y (mm)')
    end
end
sgtitle(['Reconstruction with ' num2str(maxIter) ' iterations and ' num2str(dynamic_range_dB) 'dB ' ])
