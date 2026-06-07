%% Reconstruct the images

% Envelope
V_env = abs(V_fine);

vol_true = V_fine2;

% Normalise
V_env = V_env / max(V_env(:));

% Log compression (B-mode style)
dynamic_range_dB = 15;   % typical ultrasound range
V_dB = 20*log10(V_env + eps);  % avoid log(0)

% Clip dynamic range
V_dB(V_dB < -dynamic_range_dB) = -dynamic_range_dB;

figure;

% MIP (3 projections)

eps = 1e-6;
dynamic_range_dB = 15;
vol_mag = abs(V_fine);

% MIP along z (XY view)
% For each (x,y), take the maximum value along the z direction
mip_xy1 = max(V_fine0, [], 3);
mip_xy_db1 = 20*log10(mip_xy1 / max(mip_xy1(:)) + eps);

% MIP along y (XZ view)
% For each (x,z), take the maximum value along the y direction
mip_xz1 = squeeze(max(V_fine0, [], 2));
mip_xz_db1 = 20*log10(mip_xz1 / max(mip_xz1(:)) + eps);

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
mip_xy2 = max(vol_true, [], 3);
mip_xy_db2 = 20*log10(mip_xy2 / max(mip_xy2(:)) + eps);

% MIP along y (XZ view)
% For each (x,z), take the maximum value along the y direction
mip_xz2 = squeeze(max(vol_true, [], 2));
mip_xz_db2 = 20*log10(mip_xz2 / max(mip_xz2(:)) + eps);

% plot by dimension grid
figure

t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
t.Padding = 'loose'; % gives space for labels

% =========================================================
% LEFT COLUMN — GROUND TRUTH
% =========================================================

nexttile(1)
imagesc(x_range*1e3, y_range*1e3, mip_xy_db1.')
axis tight
set(gca, 'YDir', 'reverse');
colormap gray
caxis([-dynamic_range_dB 0])
title('\bf Ground Truth')
ylabel('Elevational (mm)','FontSize',12)

nexttile(4)
imagesc(x_range*1e3, z_range*1e3, mip_xz_db1.')
axis tight
set(gca, 'YDir', 'reverse')
colormap gray
caxis([-dynamic_range_dB 0])
ylabel('Depth (mm)','FontSize',12)

nexttile(2)
imagesc(x_range*1e3, y_range*1e3, mip_xy_db2.')
axis tight
set(gca, 'YDir', 'reverse');
colormap gray
caxis([-dynamic_range_dB 0])
title('\bf Ideal Reconstruction')
%ylabel('Elevational (mm)','FontSize',12)

nexttile(5)
imagesc(x_range*1e3, z_range*1e3, mip_xz_db2.')
axis tight
set(gca, 'YDir', 'reverse')
colormap gray
caxis([-dynamic_range_dB 0])
%ylabel('Depth (mm)','FontSize',12)
% =========================================================
% RIGHT COLUMN — RECONSTRUCTED
% =========================================================

nexttile(3)
imagesc(x_range*1e3, y_range*1e3, mip_xy_db.')
axis tight
set(gca, 'YDir', 'reverse');
colormap gray
caxis([-dynamic_range_dB 0])
title('\bf Simulated Reconstruction')

nexttile(6)
imagesc(x_range*1e3, z_range*1e3, mip_xz_db.')
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
        ax(i).YLabel.Position(1) = -0.2;   % shift left (fix overlap)
        ax(i).YLabel.Position(2) = 0.5;     % vertical centre
    end
end