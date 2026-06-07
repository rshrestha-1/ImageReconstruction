function [v, labels] = carotid_phantom_tapered(Nx, Ny, Nz, scenario)

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

%% Progressive stenosis + asymmetric plaque

% Normalised z progression [0 → 1]
z_idx = linspace(0,1,Nz);
z_map = reshape(z_idx,1,1,Nz);

% 1. Circular lumen shrinkage

r_lumen_shrink = r_lumen * (1 - 0.5*z_map);  
% up to 50% narrowing

mask_blood = r <= r_lumen_shrink;

v(mask_blood) = Z_blood;
labels(mask_blood) = "blood";

% Wall updates accordingly
mask_wall = (r > r_lumen_shrink) & (r <= r_wall);
v(mask_wall) = Z_wall;
labels(mask_wall) = "wall";

% 2. Original plaque positions
lipid_center = [0.2, 0];
fibrous_center = [0.15, 0];
calc_center = [0.15, 0];

r_lipid_xy   = sqrt((x - lipid_center(1)).^2 + (y - lipid_center(2)).^2);
r_fibrous_xy = sqrt((x - fibrous_center(1)).^2 + y.^2);
r_calc_xy    = sqrt((x - calc_center(1)).^2 + y.^2);

% 3. Plaque growth (asymmetric)
fibrous_radius = 0.25 * z_map;
lipid_radius   = 0.15 * max(z_map - 0.3, 0);
calc_radius    = 0.08 * max(z_map - 0.6, 0);

fibrous_mask = (r_fibrous_xy < fibrous_radius);
lipid_mask   = (r_lipid_xy   < lipid_radius);
calc_mask    = (r_calc_xy    < calc_radius);

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
