%% 5(+). Interpolation

Nx_fine = 84;
Ny_fine = 61;
Nz_fine = 7;

x_fine = linspace(min(x_range), max(x_range), Nx_fine);
y_fine = linspace(min(y_range), max(y_range), Ny_fine);
z_fine = linspace(min(z_range), max(z_range), Nz_fine);

[xq, yq, zq] = ndgrid(x_fine, y_fine, z_fine);

V_fine0 = interpn(x_grid, y_grid, z_grid, x, xq, yq, zq, 'linear');
V_fine = interpn(x_grid, y_grid, z_grid, volume2, xq, yq, zq, 'linear');
V_fine2 = interpn(x_grid, y_grid, z_grid, vol, xq, yq, zq, 'linear');

% Envelope
V_env = abs(V_fine0);

% Normalise
V_env = V_env / max(V_env(:));

% Log compression (B-mode style)
dynamic_range_dB = 15;   % typical ultrasound range
V_dB = 20*log10(V_env + eps);  % avoid log(0)

% Clip dynamic range
V_dB(V_dB < -dynamic_range_dB) = -dynamic_range_dB;

% Plot all depth slices

figure;

for k = 1:Nz_fine
    
    subplot(ceil(Nz_fine/4), 4, k)  % 4 slices per row
    
    imagesc(x_fine*1e3, y_fine*1e3, squeeze(V_dB(:,:,k))')
    
    axis image
    set(gca,'YDir','normal')   % standard Cartesian
    colormap(gray)
    caxis([-dynamic_range_dB 0])
    
    title(sprintf('z = %.2f mm, at slice = %d', z_fine(k)*1e3, k))
    
    if k == 1
        xlabel('x (mm)')
        ylabel('y (mm)')
    end
end

sgtitle('3D LSQR Reconstruction (Log Compressed B-mode) with interpolation')
colorbar

% MIP (3 projections)

eps = 1e-6;
vol_mag = abs(V_fine);

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