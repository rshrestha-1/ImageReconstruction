%instead of repeating the bmode equations they are in one function

load('V_fine.mat');
load('V_original.mat');
load('volume5.mat');
load('x_fine.mat');
load('y_fine.mat');
load('z_range_asp.mat');

volume_bmode = @(V) ...
    max(min(20*log10(abs(V)/max(abs(V(:))) + eps), 0), -dynamic_range_dB);


volume_orig = V_original;
volume_ASA  = V_fine;
volume_PINN = volume5;

V_orig_dB = volume_bmode(volume_orig);
V_ASA_dB  = volume_bmode(volume_ASA);
V_PINN_dB = volume_bmode(volume_PINN);

% Choose 2 slices
slice_ids = [1 2];
dynamic_range_dB = 15;  


figure
t = tiledlayout(2,3,'TileSpacing','compact','Padding','compact');
t.Padding = 'loose';

titles = {'\bf Original', '\bf ASA', '\bf PINNs'};

for i = 1:2
    k = slice_ids(i);

    % =========================
    % Column 1: Original
    % =========================
    nexttile((i-1)*3 + 1)
    imagesc(x_fine*1e3, y_fine*1e3, squeeze(V_orig_dB(:,:,k))')
    axis image
    set(gca,'YDir','normal')
    colormap gray
    clim([-dynamic_range_dB 0])
    if i == 1
        title(sprintf('%s\nz = %.2f mm', titles{1}, z_range_asp(k)*1e3))
    else
        title(sprintf('z = %.2f mm', z_range_asp(k)*1e3))
    end

    if i == 1
        ylabel('Elevational (mm)','FontSize',12)
    end

    if i == 2
        ylabel('Elevational (mm)','FontSize',12)
    end
    % =========================
    % Column 2: ASA
    % =========================
    nexttile((i-1)*3 + 2)
    imagesc(x_fine*1e3, y_fine*1e3, squeeze(V_ASA_dB(:,:,k))')
    axis image
    set(gca,'YDir','normal')
    colormap gray
    clim([-dynamic_range_dB 0])

    if i == 1
        title(sprintf('%s\nz = %.2f mm', titles{2}, z_range_asp(k)*1e3))
    else
        title(sprintf('z = %.2f mm', z_range_asp(k)*1e3))
    end
    % =========================
    % Column 3: PINNs
    % =========================
    nexttile((i-1)*3 + 3)
    imagesc(x_fine*1e3, y_fine*1e3, squeeze(V_PINN_dB(:,:,k))')
    axis image
    set(gca,'YDir','normal')
    colormap gray
    clim([-dynamic_range_dB 0])
    if i == 1
        title(sprintf('%s\nz = %.2f mm', titles{3}, z_range_asp(k)*1e3))
    else
        title(sprintf('z = %.2f mm', z_range_asp(k)*1e3))
    end
end

% =========================
% SHARED LABELS
% =========================
xlabel(t,'Lateral (mm)','FontSize',13)

% =========================
% COLORBAR (GLOBAL)
% =========================
cb = colorbar;
cb.Layout.Tile = 'east';
ylabel(cb,'Intensity (dB)','FontSize',12)
cb.FontSize = 12;

% =========================
% UNIFORM AXIS FONT
% =========================
ax = findall(gcf,'Type','axes');
set(ax,'FontSize',12)

% =========================
% ALIGN Y-LABELS
% =========================
for i = 1:length(ax)
    if strcmp(ax(i).YLabel.String,'Elevational (mm)')
        ax(i).YLabel.Units = 'normalized';
        ax(i).YLabel.Position(1) = -0.25;
        ax(i).YLabel.Position(2) = 0.5;
    end
end

% =========================
% TITLE
% =========================
sgtitle(sprintf('\\bf Comparison of Original, Angular Spectrum and PINNs Reconstructions'));
