%% RESET

close all; clear all; clc

%% 1a. Set H
H = H_5;
%% 1b. Set u
u = udata_3points;
%% 4. Autoconvolution

Ne = 5;        % number of elements
Nt = 4575/Ne;     % number of time samples per element
Nv = Nx*Ny*Nz;     % number of voxels

eps = 1e-6;

Hf_3D = reshape(H, [Nt, Ne, Nv]);
H = zeros(Ne*Nt, Nv);

for v = 1:Nv
    
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
        
        H(row_start:row_end, v) = h_echo;
        
    end
    
end
%% 2. Voxel Map
% Define number of voxels
Nx = 19; % number of lateral voxels
Ny = 19; % number of height/elevational voxels
Nz = 11; % number of depth voxels

N_voxels = Nx * Ny * Nz; % Total number of voxels

% Define range of voxel map
x_range = linspace(-2.2e-3, 2.2e-3, Nx); % mm
y_range = linspace(-2.2e-3, 2.2e-3, Ny); % mm 
z_range = linspace(11.24e-3, 12.44e-3, Nz); % mm

[x_grid, y_grid, z_grid] = ndgrid(x_range, y_range, z_range);

voxel_positions = [x_grid(:), y_grid(:), z_grid(:)];
%% 3a. Normalise H
H = H / max(abs(H(:)));
%% 3b. Normalise u
u = u / max(abs(u));
%% 5. LSQR
tol = 1e-6;
maxit = 23;

v = lsqr(Hn,u,tol,maxit);

volume = reshape(v,Nx,Ny,Nz);
%% Scaled volume viewer
vmax = max(volume(:))*0.2;

volume_scaled = volume;
volume_scaled(volume_scaled > vmax) = vmax;  % clip high values

volumeViewer(volume_scaled)
%% Tikhonov Regularisation version

% Add small regularization to stabilize
lambda = 1e-3;
[v,flag] = lsqr([H; lambda*eye(size(H,2))],[u; zeros(size(H,2),1)],1e-4,23);

volume = reshape(v,Nx,Ny,Nz);
%% 6. Volume show (percentile)

figure
vmin = prctile(volume(:),1);
vmax = prctile(volume(:),99);

montage(volume,'DisplayRange',[vmin vmax])
title('Reconstructed slices')
colormap gray
%% Lower percentage of Dynamic Range
figure
montage(volume,'DisplayRange',[0 max(volume(:))*0.2])
%% 3D B-mode style visualisation (all XY slices)

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

% MIP along z (XY view)
% For each (x,y), take the maximum value along the z direction
mip_xy = max(volume, [], 3);

% MIP along y (XZ view)
% For each (x,z), take the maximum value along the y direction
mip_xz = squeeze(max(volume, [], 2));

% MIP along x (YZ view)
% For each (y,z), take the maximum value along the x direction
mip_yz = squeeze(max(volume, [], 1));

% plot
figure

% XY projection (looking down the z-axis)
subplot(1,3,1)
imagesc(mip_xy)
axis image
colormap gray
colorbar
title('MIP - XY view (projection along Z)')


% XZ projection (looking along the y-axis)
subplot(1,3,2)
imagesc(mip_xz)
axis image
colormap gray
colorbar
title('MIP - XZ view (projection along Y)')


% YZ projection (looking along the x-axis)
subplot(1,3,3)
imagesc(mip_yz)
axis image
colormap gray
colorbar
title('MIP - YZ view (projection along X)')

sgtitle('Maximum Intensity Projections of 3D Reconstruction')
%% Visualise voxel sensitity
sens = reshape(sum(abs(H),1),Nx,Ny,Nz);
figure;
montage(sens,'DisplayRange',[])
%% Test average column
col_norms = sqrt(sum(H.^2,1));
vol = reshape(col_norms,Nx,Ny,Nz);
montage(vol,'DisplayRange',[])
%% Column Normalisation
col_norms = sqrt(sum(H.^2,1));
Hn = H ./ col_norms;

sensn = reshape(sum(abs(Hn),1),Nx,Ny,Nz);

figure;
montage(sensn,'DisplayRange',[])
title('Voxel sensitivity after column normalisation')
%% Test column correlation

rand_vox = randperm(Nx*Ny*Nz,200);

C = corrcoef(H(:,rand_vox));
figure;
imagesc(C)
colorbar
%% Perfectly simulated u

x = zeros(Nx,Ny,Nz);

% If Voxel Number known
%[i, j, k] = ind2sub([Nx, Ny, Nz], 2708);

% If slices are known
idx = sub2ind([Nx, Ny, Nz], 20, 9, 8);
x(idx) = 1;

% Else slices known
x(10,5,10) = 1;
x(10,10,10) = 1;
x(10,15,10) = 1;

b = H * x(:);

xrec = lsqr(H,b,1e-6,15);

figure;
vol = reshape(xrec,Nx,Ny,Nz);
montage(vol,'DisplayRange',[0 max(vol(:))*0.2])

%% 3D B-mode style visualisation (all XY slices)

% Envelope
V_env = abs(vol);

% Normalise
V_env = V_env / max(V_env(:));

% Log compression (B-mode style)
dynamic_range_dB = 15;   % 0 to 60 typical ultrasound range
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
dynamic_range_dB = 18;
vol_mag = abs(vol);

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

%% Plot H (autoconvolved)

N_voxel = 3335;
Nt = 4575/5;

signal_voxel1 = H(1:Nt, N_voxel);   
signal_voxel2 = H(Nt+1:2*Nt, N_voxel); 
signal_voxel3 = H(2*Nt+1:3*Nt, N_voxel); 
t = 1:Nt;     

figure;

% Element 1
subplot(1,3,1);
plot(t, signal_voxel1, 'LineWidth', 1.5);
xlabel('Time Samples');
ylabel('Amplitude');
title(['H Signal for Voxel ', num2str(N_voxel)],'with element 1');
grid on;

% Element 2
subplot(1,3,2);
plot(t, signal_voxel2, 'LineWidth', 1.5);
xlabel('Time Samples');
ylabel('Amplitude');
title(['H Signal for Voxel ', num2str(N_voxel)],'with element 2');
grid on;

% Element 3
subplot(1,3,3);
plot(t, signal_voxel3, 'LineWidth', 1.5);
xlabel('Time Samples');
ylabel('Amplitude');
title(['H Signal for Voxel ', num2str(N_voxel)],'with element 3');
grid on;
%% Plot u

U = udata_3points;

signal_u1 = U(1:Nt);   
signal_u2 = U(Nt+1:2*Nt); 
signal_u3 = U(2*Nt+1:3*Nt); 

t = 1:Nt;      

figure;

% Element 1
subplot(1,3,1);
plot(t, signal_u1, 'LineWidth', 1.5);
xlabel('Time Samples');
ylabel('Amplitude');
title('u Signal with element 1');
grid on;

% Element 2
subplot(1,3,2);
plot(t, signal_u2, 'LineWidth', 1.5);
xlabel('Time Samples');
ylabel('Amplitude');
title('u Signal with element 2');
grid on;

% Element 3
subplot(1,3,3);
plot(t, signal_u3, 'LineWidth', 1.5);
xlabel('Time Samples');
ylabel('Amplitude');
title('u Signal with element 3');
grid on;
%% Plot v
figure;
plot(v, 'LineWidth', 1.5);
xlabel('Index');
ylabel('Intensity');
title('Intensity Values');
grid on;
