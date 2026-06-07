%% Reconstruct the images

% Original
% Envelope
V_env = abs(volume);

% Normalise
V_env = V_env / max(V_env(:));

% Log compression (B-mode style)
dynamic_range_dB = 15;   % typical ultrasound range
V_dB = 20*log10(V_env + eps);  % avoid log(0)

% Clip dynamic range
V_dB(V_dB < -dynamic_range_dB) = -dynamic_range_dB;

% Downsampled 2
V_2 = abs(volume2);

% Downsampled 2
V_3 = abs(volume3);

% Downsampled 2
V_4 = abs(volume4);

% Downsampled 2
V_5 = abs(volume5);

% MIP (3 projections)

eps = 1e-6;
dynamic_range_dB = 10;
vol_mag = abs(volume);

% Original
% MIP along z (XY view)
% For each (x,y), take the maximum value along the z direction
mip_xy = max(vol_mag, [], 3);
mip_xy_db = 20*log10(mip_xy / max(mip_xy(:)) + eps);

% MIP along y (XZ view)
% For each (x,z), take the maximum value along the y direction
mip_xz = squeeze(max(vol_mag, [], 2));
mip_xz_db = 20*log10(mip_xz / max(mip_xz(:)) + eps);

% Factor 2
% MIP along z (XY view)
% For each (x,y), take the maximum value along the z direction
mip_xy2 = max(V_2, [], 3);
mip_xy_db2 = 20*log10(mip_xy2 / max(mip_xy2(:)) + eps);

% MIP along y (XZ view)
% For each (x,z), take the maximum value along the y direction
mip_xz2 = squeeze(max(V_2, [], 2));
mip_xz_db2 = 20*log10(mip_xz2 / max(mip_xz2(:)) + eps);

% Factor 3
% MIP along z (XY view)
% For each (x,y), take the maximum value along the z direction
mip_xy3 = max(V_3, [], 3);
mip_xy_db3 = 20*log10(mip_xy3 / max(mip_xy3(:)) + eps);

% MIP along y (XZ view)
% For each (x,z), take the maximum value along the y direction
mip_xz3 = squeeze(max(V_3, [], 2));
mip_xz_db3 = 20*log10(mip_xz3 / max(mip_xz3(:)) + eps);

% Factor 4
% MIP along z (XY view)
% For each (x,y), take the maximum value along the z direction
mip_xy4 = max(V_4, [], 3);
mip_xy_db4 = 20*log10(mip_xy4 / max(mip_xy4(:)) + eps);

% MIP along y (XZ view)
% For each (x,z), take the maximum value along the y direction
mip_xz4 = squeeze(max(V_4, [], 2));
mip_xz_db4 = 20*log10(mip_xz4 / max(mip_xz4(:)) + eps);

% Factor 5
% MIP along z (XY view)
% For each (x,y), take the maximum value along the z direction
mip_xy5 = max(V_5, [], 3);
mip_xy_db5 = 20*log10(mip_xy5 / max(mip_xy5(:)) + eps);

% MIP along y (XZ view)
% For each (x,z), take the maximum value along the y direction
mip_xz5 = squeeze(max(V_5, [], 2));
mip_xz_db5 = 20*log10(mip_xz5 / max(mip_xz5(:)) + eps);

% plot by dimension grid
figure

t = tiledlayout(2,5,'TileSpacing','compact','Padding','compact');
t.Padding = 'loose'; % gives space for labels

% Original

nexttile(1)
imagesc(x_range*1e3, y_range*1e3, mip_xy_db.')
axis tight
set(gca, 'YDir', 'normal');
colormap gray
caxis([-dynamic_range_dB 0])
title('\bf Original')
ylabel('Elevational (mm)','FontSize',12)

nexttile(6)
imagesc(x_range*1e3, z_range*1e3, mip_xz_db.')
axis tight
set(gca, 'YDir', 'reverse')
colormap gray
caxis([-dynamic_range_dB 0])
ylabel('Depth (mm)','FontSize',12)

% Downsample 2

nexttile(2)
imagesc(x_range*1e3, y_range*1e3, mip_xy_db2.')
axis tight
set(gca, 'YDir', 'normal');
colormap gray
caxis([-dynamic_range_dB 0])
title('\bf DS ×2')

nexttile(7)
imagesc(x_range*1e3, z_range*1e3, mip_xz_db2.')
axis tight
set(gca, 'YDir', 'reverse')
colormap gray
caxis([-dynamic_range_dB 0])

% Downsample 3

nexttile(3)
imagesc(x_range*1e3, y_range*1e3, mip_xy_db3.')
axis tight
set(gca, 'YDir', 'normal');
colormap gray
caxis([-dynamic_range_dB 0])
title('\bf DS ×3')

nexttile(8)
imagesc(x_range*1e3, z_range*1e3, mip_xz_db3.')
axis tight
set(gca, 'YDir', 'reverse')
colormap gray
caxis([-dynamic_range_dB 0])

% Downsample 4

nexttile(4)
imagesc(x_range*1e3, y_range*1e3, mip_xy_db4.')
axis tight
set(gca, 'YDir', 'normal');
colormap gray
caxis([-dynamic_range_dB 0])
title('\bf DS ×4')

nexttile(9)
imagesc(x_range*1e3, z_range*1e3, mip_xz_db4.')
axis tight
set(gca, 'YDir', 'reverse')
colormap gray
caxis([-dynamic_range_dB 0])

% Downsample 5

nexttile(5)
imagesc(x_range*1e3, y_range*1e3, mip_xy_db5.')
axis tight
set(gca, 'YDir', 'normal');
colormap gray
caxis([-dynamic_range_dB 0])
title('\bf DS ×5')

nexttile(10)
imagesc(x_range*1e3, z_range*1e3, mip_xz_db5.')
axis tight
set(gca, 'YDir', 'reverse')
colormap gray
caxis([-dynamic_range_dB 0])

% =========================================================
% SHARED LABELS + FORMATTING
% =========================================================

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

% =========================================================
% FINAL VISUAL TIDYING (ALIGN Y LABELS)
% =========================================================
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

sgtitle('\bf Comparison of Original and Downsampled–Interpolated Reconstructions', ...
    'Interpreter','tex');