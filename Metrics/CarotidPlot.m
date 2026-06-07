%% Reconstruct the images

% Define range of voxel map
z_range = linspace(-10e-3, 10e-3, Nx); % m
x_range = linspace(-16e-3, 16e-3, Ny); % m
y_range = linspace(8e-3, 40e-3, Nz); % m

% Envelope
V_env = abs(v_true_3D);

% Normalise
V_env = V_env / max(V_env(:));

% Log compression (B-mode style)
dynamic_range_dB = 15;   % typical ultrasound range
V_dB = 20*log10(V_env + eps);  % avoid log(0)

% Clip dynamic range
V_dB(V_dB < -dynamic_range_dB) = -dynamic_range_dB;

% MIP (3 projections)

eps = 1e-6;
dynamic_range_dB = 15;
voltrueside = abs(v_true_3D2);
volrecxy = abs(v_rec_3D);
volreczy = abs(v_rec_3D2);

% MIP zy view
mip_yz = squeeze(max(voltrueside, [], 1));
mip_yz_db = 20*log10(mip_yz / max(mip_yz(:)) + eps);

% MIP xy view Recon
mip_xy1 = max(volrecxy, [], 3);
mip_xy_db1 = 20*log10(mip_xy1 / max(mip_xy1(:)) + eps);

% MIP zy view Recon
mip_yz2 = squeeze(max(volreczy, [], 1));
mip_yz_db2 = 20*log10(mip_yz2 / max(mip_yz2(:)) + eps);
 
figure;

t = tiledlayout(2,2,'TileSpacing','compact','Padding','loose');

% =========================================================
% TOP ROW — GROUND TRUTH
% =========================================================

% --- YZ (left) ---
% --- First slice (z = 1) instead of YZ MIP ---
ax1 = nexttile(1);
slice_xy = squeeze(V_dB(:,:,1));   % first depth slice
imagesc(x_range*1e3, y_range*1e3, slice_xy.')

axis tight
set(gca,'YDir','reverse')
colormap gray
caxis([-dynamic_range_dB 0])

ylabel('Depth (mm)','FontSize',12)
% no xlabel (top row)

% --- XZ (right) ---
ax2 = nexttile(2);
imagesc(z_range*1e3, y_range*1e3, mip_yz_db)
axis tight
set(gca,'YDir','reverse')
colormap gray
caxis([-dynamic_range_dB 0])

% no labels here (aligned with left)

% Row title
%text(0.5,1.05,'\bf Ground Truth','Units','normalized', ...
    %'HorizontalAlignment','center','Parent',ax1,'FontSize',14)

% =========================================================
% BOTTOM ROW — RECONSTRUCTION
% =========================================================

% --- YZ (left) ---
ax3 = nexttile(3);
imagesc(x_range*1e3, y_range*1e3, mip_xy_db1.')

axis tight
set(gca,'YDir','reverse')
colormap gray
caxis([-dynamic_range_dB 0])

xlabel('Elevational (mm)','FontSize',12)
ylabel('Depth (mm)','FontSize',12)

% --- XZ (right) ---
ax4 = nexttile(4);
imagesc(z_range*1e3, y_range*1e3, mip_yz_db2)
axis tight
set(gca,'YDir','reverse')
colormap gray
caxis([-dynamic_range_dB 0])

xlabel('Lateral (mm)','FontSize',12)

% no ylabel

% Row title
%text(0.5,1.05,'\bf Reconstruction','Units','normalized', ...
    %'HorizontalAlignment','center','Parent',ax3,'FontSize',14)

% =========================================================
% COLORBAR (FULL HEIGHT)
% =========================================================

cb = colorbar;
cb.Layout.Tile = 'east';
ylabel(cb,'Intensity (dB)','FontSize',12)
cb.FontSize = 12;

% =========================================================
% ALIGN Y-LABELS NICELY (LEFT COLUMN ONLY)
% =========================================================

for ax = [ax1, ax3]
    ax.YLabel.Units = 'normalized';
    ax.YLabel.Position(1) = -0.18;
    ax.YLabel.Position(2) = 0.5;
end

% =========================================================
% FONT CONSISTENCY
% =========================================================

all_axes = findall(gcf,'Type','axes');
set(all_axes,'FontSize',12)

% =========================================================
% ROW TITLES (FIXED & STABLE)
% =========================================================

% Get positions (normalized to figure)
pos1 = ax1.Position; % top-left
pos3 = ax3.Position; % bottom-left

% Top of top row
top_row_top = pos1(2) + pos1(4);

% Top of bottom row
bottom_row_top = pos3(2) + pos3(4);

drawnow;

pos1 = ax1.Position;
pos3 = ax3.Position;

top_row_top = pos1(2) + pos1(4);
bottom_row_top = pos3(2) + pos3(4);

gap_center = (top_row_top + bottom_row_top) / 2;

% =========================================================
% TITLES (TINY CONTROLLED GAP)
% =========================================================

% Ground Truth (slightly above top row)
annotation('textbox', ...
    [0.5-0.15, top_row_top, 0.3, 0.05], ...
    'String','\bf Ground Truth Phantom', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','center', ...
    'FontSize',14);

annotation('textbox', ...
    [0.5-0.25, gap_center - 0.225, 0.5, 0.05], ...
    'String','\bf Reconstructed Phantom', ...
    'EdgeColor','none', ...
    'HorizontalAlignment','center', ...
    'FontSize',14, ...
    'FitBoxToText','off');
% =========================================================
drawnow;

% =========================================================
% CREATE SMALL EXTRA ROW GAP (SAFE METHOD)
% =========================================================

shrink = 0.92;  % <-- control gap here (0.90–0.95 recommended)

for ax = [ax1, ax2, ax3, ax4]
    p = ax.Position;
    
    % shrink height slightly
    new_h = p(4) * shrink;
    
    % keep bottom anchored
    p(2) = p(2) + (p(4) - new_h);
    p(4) = new_h;
    
    ax.Position = p;
end
% =========================================================
% ADJUST COLUMN WIDTHS (LEFT SMALLER, RIGHT LARGER)
% =========================================================

scale_left  = 0.6;  % shrink left plots (try 0.8–0.9)
scale_right = 1.50;  % expand right plots (try 1.05–1.15)

% --- LEFT COLUMN (ax1, ax3) ---
for ax = [ax1, ax3]
    pos = ax.Position;
    
    % shrink width
    new_width = pos(3) * scale_left;
    
    % keep left edge fixed
    pos(3) = new_width;
    
    ax.Position = pos;
end

% --- RIGHT COLUMN (ax2, ax4) ---
for ax = [ax2, ax4]
    pos = ax.Position;
    
    % expand width
    delta = pos(3) * (scale_right - 1);
    
    % shift left slightly so it expands inward
    pos(1) = pos(1) - delta;
    pos(3) = pos(3) * scale_right;
    
    ax.Position = pos;
end
