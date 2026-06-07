%% Perfectly simulated u

x = zeros(Nx,Ny,Nz);

% If Voxel Number known
%[i, j, k] = ind2sub([Nx,10 Ny, Nz], 2708);

% If slices are known
%idx = sub2ind([Nx, Ny, Nz], 20, 9, 8);
%x(idx) = 1;

% Else slices known
x(10,5,4) = 1;
x(10,10,4) = 1;
x(10,15,4) = 1;
%% 

b = H * x(:);

xrec = lsqr(Hn,b,1e-6,30);

% Threshold
xrec(abs(xrec) < 0.1*max(abs(xrec))) = 0;

vol = reshape(xrec,Nx,Ny,Nz);

%% H_30 u

v = zeros(Nx, Ny, Nz);

x_range = 24:31;
z_range = 2:5;

y_ranges = [10:13, 30:33, 50:53];

v(x_range, y_ranges, z_range) = 1;
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
    
    %title(sprintf('z = %.2f mm, at slice = %d', z_range(k)*1e3, k))
    
    if k == 1
        xlabel('x (mm)')
        ylabel('y (mm)')
    end
end

sgtitle('3D LSQR Reconstruction (Log Compressed B-mode)')
colorbar
%% MIP (3 projections)

eps = 1e-6;
dynamic_range_dB = 15;
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

U = b;

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
plot(xrec, 'LineWidth', 1.5);
xlabel('Index');
ylabel('Intensity');
title('Intensity Values');
grid on;