function [v, labels] = carotid_phantom(Nx, Ny, Nz, scenario)

[x, y, z] = ndgrid(linspace(-1,1,Nx), linspace(-1,1,Ny), linspace(-1,1,Nz));

v = zeros(Nx,Ny,Nz);
labels = strings(Nx,Ny,Nz);

%% Acoustic impedance
Z_blood   = 1060 * 1570;
Z_wall    = 1100 * 1600;
Z_lipid   = 920  * 1450;
Z_fibrous = 1150 * 1670;
Z_calc    = 1500 * 3500;

%% Geometry
r_lumen = 0.4;
r_wall  = 0.5;
r = sqrt(x.^2 + y.^2);

%% Base anatomy
mask_blood = r <= r_lumen;
v(mask_blood) = Z_blood;
labels(mask_blood) = "blood";

mask_wall = (r > r_lumen) & (r <= r_wall);
v(mask_wall) = Z_wall;
labels(mask_wall) = "wall";

%% Plaque geometry (same shapes reused)
lipid_center = [0.2, 0, 0];
lipid_mask = sqrt((x-lipid_center(1)).^2 + (y-lipid_center(2)).^2 + z.^2) < 0.15;

fibrous_mask = sqrt((x-0.15).^2 + y.^2 + z.^2) < 0.2;

% calc_mask = sqrt((x+0.2).^2 + y.^2 + z.^2) < 0.1;
calc_mask = sqrt((x).^2 + y.^2 + z.^2) < 0.1;

%% Scenario switch
switch lower(scenario)

    case "healthy"
        % nothing extra

    case "lipid"
        v(lipid_mask) = Z_lipid;
        labels(lipid_mask) = "lipid";

    case "fibrous"
        v(fibrous_mask) = Z_fibrous;
        labels(fibrous_mask) = "fibrous";

    case "calcified"
        v(calc_mask) = Z_calc;
        labels(calc_mask) = "calcified";

    case "mixed"
        v(lipid_mask) = Z_lipid;
        labels(lipid_mask) = "lipid";

        v(fibrous_mask & ~lipid_mask) = Z_fibrous;
        labels(fibrous_mask & ~lipid_mask) = "fibrous";

        v(calc_mask) = Z_calc;
        labels(calc_mask) = "calcified";

    otherwise
        error("Unknown scenario type");

end

end