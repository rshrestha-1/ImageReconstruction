%% Sanity check: matched-filter 3D imaging visualisation

clear; clc;
%% Load
% point scatterer loads
load("/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/H_maskless.mat")
load("/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/u_maskless.mat")
%% Load
% point scatterer loads
load("/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/H_maskless2.mat")
load("/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/u_maskless2.mat")
load("/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/u_maskless3.mat")
load("/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/u_maskless4.mat")
load("/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/I_new.mat")
%% Load
% old h and logo
load("/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/udata_logo.mat")
load("/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/H_f.mat")
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/udata_point.mat')
%% Load
% new correct maskless
load("/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/H_maskless3.mat")
load("/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/correct_u_maskless.mat")
load("/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/diagonal_u_maskless.mat")
load("/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/I_maskless.mat")
load("/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/I_point_u_maskless.mat")
load("/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/maskless_u_1790.mat")
%% Load
% Experimental 1
load("/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/H_exp.mat")
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/u_exp.mat')

%% Load
% Mask sim
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/u_data_slice8.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/H_293.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/H_303.mat')
%% Load
% Experimental 2
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/u_exp_new.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/H_experimental_new.mat')
%% Load
% 5 element H
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/udata_3points.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/H_5.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/trial_1.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/trial2.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/trial3.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/u_data_diagonal.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/timegated_udata_d.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/udata_3_middle.mat')
%% Load
% 15 element H
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/u_3_15elem.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/H_15.mat')

%% Load
% 15 element H

%% Voxel Map
% Define number of voxels
Nx = 19; % number of lateral voxels
Ny = 19; % number of height/elevational voxels
Nz = 8; % number of depth voxels

N_voxels = Nx * Ny * Nz; % Total number of voxels

% Define range of voxel map
x_range = linspace(-2.18e-3, 2.18e-3, Nx); % mm
y_range = linspace(-2.18e-3, 2.18e-3, Ny); % mm 
z_range = linspace(5.04e-3, 6.69e-3, Nz); % mm

[x_grid, y_grid, z_grid] = ndgrid(x_range, y_range, z_range);

voxel_positions = [x_grid(:), y_grid(:), z_grid(:)];

%% Simulation parameters
N_elements = 3; % Number of elements in the transducer
N_t = 1830; % Number of time samples for the measurement

%% Simulate calibration matrix H
% Each column = flattened response of a voxel (N_elements*N_t x 1)
H = randn(N_elements*N_t, N_voxels) * 0.1;  % small random noise

% Autoconvolve H externally

%% Simulate measurement u using the true voxel only
true_voxel_index = 4846;

v_true = zeros(N_voxels,1);
v_true(true_voxel_index) = 1;   % strong unit reflector

% Make that voxel column stronger
H(:,true_voxel_index) = randn(N_elements*N_t,1);

u = H * v_true;
%% Reconstruct with LSQR

max_iter = 100;
tol = 1e-8;

v = lsqr(H, u, tol, max_iter);

true_x = voxel_positions(true_voxel_index,1);
true_y = voxel_positions(true_voxel_index,2);
true_z = voxel_positions(true_voxel_index,3);

fprintf('True voxel %d at x = %.4f mm, y = %.4f mm, z = %.4f mm\n', ...
    true_voxel_index, true_x*1e3, true_y*1e3, true_z*1e3);
%% Check peak location

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

% Envelope
V_env = abs(V);

% Normalise
V_env = V_env / max(V_env(:));

% Log compression (B-mode style)
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


%% Relative reconstruction error
% fprintf('\nTrue index: %d\n', true_voxel_index);
fprintf('Recovered index: %d\n', recovered_index);

error_norm = norm(v - v_true) / norm(v_true);
fprintf('Relative reconstruction error: %.4e\n', error_norm);

%% REAL

%% Voxel Map
% Define number of voxels
Nx = 19; % number of lateral voxels
Ny = 19; % number of height/elevational voxels
Nz = 8; % number of depth voxels

N_voxels = Nx * Ny * Nz; % Total number of voxels

% Define range of voxel map
x_range = linspace(-2.18e-3, 2.18e-3, Nx); % mm
y_range = linspace(-2.18e-3, 2.18e-3, Ny); % mm 
%z_range = linspace(7.92e-3, 10.08e-3, Nz); % mm
z_range = linspace(5.04e-3, 6.69e-3, Nz); % mm

[x_grid, y_grid, z_grid] = ndgrid(x_range, y_range, z_range);

voxel_positions = [x_grid(:), y_grid(:), z_grid(:)];

%% Non-simulated calibration matrix H

H = H_f;

% Autoconvolve H externally

%% Non-simulated measurement u
u = u_maskless;
%% Reconstruct with LSQR

max_iter = 20;
tol = 1e-6;

v = lsqr(H, u, tol, max_iter);
%% Check true voxel
true_voxel_index = 3069;

true_x = voxel_positions(true_voxel_index,1);
true_y = voxel_positions(true_voxel_index,2);
true_z = voxel_positions(true_voxel_index,3);

fprintf('True voxel %d at x = %.4f mm, y = %.4f mm, z = %.4f mm\n', ...
    true_voxel_index, true_x*1e3, true_y*1e3, true_z*1e3);
%% Check peak location

[peak_val, recovered_index] = max(abs(v));

fprintf('True voxel index:      %d\n', true_voxel_index);
fprintf('Recovered voxel index: %d\n', recovered_index);

rec_x = voxel_positions(recovered_index,1);
rec_y = voxel_positions(recovered_index,2);
rec_z = voxel_positions(recovered_index,3);

fprintf('Recovered at x=%.3f mm, y=%.3f mm, z=%.3f mm\n', ...
        rec_x*1e3, rec_y*1e3, rec_z*1e3);

sorted_vals = sort(abs(v), 'descend');
fprintf('Peak-to-second ratio: %.2f\n', sorted_vals(1)/sorted_vals(2));

if recovered_index == true_voxel_index
    disp('SUCCESS: Correct voxel recovered.');
else
    disp('WARNING: Incorrect voxel recovered.');
end
%% 

V = reshape(abs(v), Nx, Ny, Nz);

%% 3D B-mode style visualisation (all XY slices)

% Envelope
V_env = abs(V);

% Normalise
V_env = V_env / max(V_env(:));

% Log compression (B-mode style)
dynamic_range_dB = 20;   % typical ultrasound range
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

%% 3D Layered XY Slices (Stacked B-mode)

% Envelope
V_env = abs(V);

% Normalise
V_env = V_env / max(V_env(:));

% Log compression
dynamic_range_dB = 9;
V_dB = 20*log10(V_env + eps);
V_dB(V_dB < -dynamic_range_dB) = -dynamic_range_dB;

figure;
hold on;

colormap(gray);
caxis([-dynamic_range_dB 0]);

for k = 1:Nz
    
    % Create constant Z plane
    Z_plane = ones(Nx, Ny) * z_range(k)*1e3;   % in mm
    
    % Plot slice as surface
    surf(x_range*1e3, ...
         y_range*1e3, ...
         Z_plane', ...                     % transpose for orientation
         squeeze(V_dB(:,:,k))', ...
         'EdgeColor','none');
end

xlabel('x (mm)')
ylabel('y (mm)')
zlabel('z (mm)')
title('3D LSQR Reconstruction (Layered B-mode)')

colorbar
view(3)
axis tight
grid on
