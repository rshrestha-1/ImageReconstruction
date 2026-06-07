%% Sanity check: matched-filter 3D imaging visualisation

clear; clc;

%% 1. Simulation parameters
% Define number of voxels
Nx = 38; % number of lateral voxels
Ny = 17; % number of height/elevational voxels
Nz = 8; % number of depth voxels

N_voxels = Nx * Ny * Nz; % Total number of voxels

% Define range of voxel map
x_range = linspace(-4.5e-3, 4.5e-3, Nx); % mm
y_range = linspace(-2e-3, 2e-3, Ny); % mm 
z_range = linspace(7.04e-3, 8.8e-3, Nz); % mm

[x_grid, y_grid, z_grid] = ndgrid(x_range, y_range, z_range);

voxel_positions = [x_grid(:), y_grid(:), z_grid(:)];

%% Simulation parameters
N_elements = 3; % Number of elements in the transducer
N_t = 1830; % Number of time samples for the measurement

%% 4. Simulate calibration matrix H
% Each column = flattened response of a voxel (N_elements*N_t x 1)
H = randn(N_elements*N_t, N_voxels) * 0.1;  % small random noise

% Autoconvolve H
H = ifft2( fft2(H) .* fft2(H) );
H = real(H);
%% 5. Simulate measurement u using the true voxel only
true_voxel_index = 513;

u_true = zeros(N_voxels,1);
u_true(true_voxel_index) = 1;   % strong unit reflector

% Make that voxel column stronger
H(:,true_voxel_index) = randn(N_elements*N_t,1);

u = H * u_true;
%% 5. Reconstruct with LSQR

max_iter = 100;
tol = 1e-8;

v = lsqr(H, u, tol, max_iter);

true_x = voxel_positions(true_voxel_index,1);
true_y = voxel_positions(true_voxel_index,2);
true_z = voxel_positions(true_voxel_index,3);

fprintf('True voxel %d at x = %.4f mm, y = %.4f mm, z = %.4f mm\n', ...
    true_voxel_index, true_x*1e3, true_y*1e3, true_z*1e3);
%% 6. Check peak location

[~, recovered_index] = max(abs(v));

fprintf('True voxel index:      %d\n', true_voxel_index);
fprintf('Recovered voxel index: %d\n', recovered_index);

if recovered_index == true_voxel_index
    disp('SUCCESS: Correct voxel recovered.');
else
    disp('WARNING: Incorrect voxel recovered.');
end

V = reshape(abs(v), Nx, Ny, Nz);

%% 3D B-mode style visualisation (all XY slices)

% --- Envelope (magnitude already used, but ensure positive)
V_env = abs(V);

% --- Normalise
V_env = V_env / max(V_env(:));

% --- Log compression (B-mode style)
dynamic_range_dB = 40;   % typical ultrasound range
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
%% VolumeViwer

volumeViewer(V_dB)


%% 
% fprintf('\nTrue index: %d\n', true_voxel_index);
fprintf('Recovered index: %d\n', recovered_index);

error_norm = norm(v - u_true) / norm(u_true);
fprintf('Relative reconstruction error: %.4e\n', error_norm);
