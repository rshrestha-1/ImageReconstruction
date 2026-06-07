%% Reconstruct the images

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

% MIP (3 projections)

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

% MIP along z (XY view)
% For each (x,y), take the maximum value along the z direction
mip_xy2 = max(vol, [], 3);
mip_xy_db2 = 20*log10(mip_xy / max(mip_xy(:)) + eps);

% MIP along y (XZ view)
% For each (x,z), take the maximum value along the y direction
mip_xz2 = squeeze(max(vol, [], 2));
mip_xz_db2 = 20*log10(mip_xz / max(mip_xz(:)) + eps);

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

t = tiledlayout(2,2,'TileSpacing','compact','Padding','compact');
t.Padding = 'loose'; % gives space for labels

% LEFT COLUMN — GROUND TRUTH

nexttile(1)
imagesc(x_range*1e3, y_range*1e3, mip_xy_db.')
axis tight
set(gca, 'YDir', 'normal');
colormap gray
caxis([-dynamic_range_dB 0])
title('\bf Ground Truth')
ylabel('Elevational (mm)','FontSize',12)

nexttile(3)
imagesc(x_range*1e3, z_range*1e3, mip_xz_db.')
axis tight
set(gca, 'YDir', 'reverse')
colormap gray
caxis([-dynamic_range_dB 0])
ylabel('Depth (mm)','FontSize',12)

% RIGHT COLUMN — RECONSTRUCTED

nexttile(2)
imagesc(x_range*1e3, y_range*1e3, mip_xy_db2.')
axis tight
set(gca, 'YDir', 'normal');
colormap gray
caxis([-dynamic_range_dB 0])
title('\bf Reconstructed')

nexttile(4)
imagesc(x_range*1e3, z_range*1e3, mip_xz_db2.')
axis tight
set(gca, 'YDir', 'reverse')
colormap gray
caxis([-dynamic_range_dB 0])

% SHARED LABELS + FORMATTING

% Shared X label
xlabel(t,'Lateral (mm)','FontSize',13)

% Colorbar (full height, right side)
cb = colorbar;
cb.Layout.Tile = 'east';
ylabel(cb,'Intensity (dB)','FontSize',12)
cb.FontSize = 12;

% Uniform font size for axes
ax = findall(gcf,'Type','axes');
set(ax,'FontSize',12)

ax = findall(gcf,'Type','axes');

for i = 1:length(ax)
    
    if strcmp(ax(i).YLabel.String,'Depth (mm)') || ...
       strcmp(ax(i).YLabel.String,'Elevational (mm)')
        
        ax(i).YLabel.Units = 'normalized';
        
        % consistent clean placement
        ax(i).YLabel.Position(1) = -0.25;   % shift left (fix overlap)
        ax(i).YLabel.Position(2) = 0.5;     % vertical centre
    end
end
