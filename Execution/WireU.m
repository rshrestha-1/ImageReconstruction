%% FIRST HAVE u and H in workspace
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/H_slice.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/H_f.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/H_fresh.mat')

load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/u_data_filtered.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/udata_logo.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/udata_point.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/udata_point2.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/u_slice21.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/u_slice18.mat')
%% 
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/H_maskless.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/H_maskless2.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/u_maskless.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/u_maskless2.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/u_maskless3.mat')
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/u_maskless4.mat')

%% Solution via LSQR

maxIter = 18;
tol = 1e-6;
u = udata_logo;

v = lsqr(H,u,tol,maxIter); % lsqr solution at 18 iterations
%% Quick autoconvolution

H_conv = ifft2( fft2(H_f) .* fft2(H_f) );
H_conv = real(H_conv);

v = lsqr(H_conv, udata_logo, tol, maxIter);
%% Filter image u

fs = 30e6; % sampling frequency        
f_l  = 3e6;      
f_h = 12e6;      

bpFilt = designfilt('bandpassiir', 'FilterOrder', 6, 'HalfPowerFrequency1', f_l, 'HalfPowerFrequency2', f_h, 'SampleRate', fs);

u_filt = filtfilt(bpFilt, u); % apply bandpass filter
%% Define voxel coordinates

% Define number of voxels
Nx = 38; % number of lateral voxels
Ny = 17; % number of height/elevational voxels
Nz = 8; % number of depth voxels

N_voxels = Nx * Ny * Nz; % Total number of voxels

% Define range of voxel map
x_range = linspace(-4.5e-3, 4.5e-3, Nx); % mm
y_range = linspace(-2e-3, 2e-3, Ny); % mm 
z_range = linspace(7.04e-3, 8.72e-3, Nz); % mm

[x_grid, y_grid, z_grid] = meshgrid(x_range, y_range, z_range); % construct voxel map

% Permute so x → y → z
x_grid = permute(x_grid, [2 1 3]);
y_grid = permute(y_grid, [2 1 3]);
z_grid = permute(z_grid, [2 1 3]);

% Flatten to voxel positions
voxel_positions = [x_grid(:), y_grid(:), z_grid(:)];

%% 2D slices layered

% Reshape v image voxel data
v_vol = permute(reshape(v, Ny, Nx, Nz), [2 1 3]); % preserve order x,y,z

% Filter per depth value for each x,y voxel
for ix = 1:Nx
    for iy = 1:Ny
        % Extract depth trace at this voxel (x,y)
        trace = squeeze(v_vol(ix, iy, :));  
        % Apply zero-phase bandpass filter
        %v_vol(ix, iy, :) = filtfilt(bpFilt, trace);  
    end
end

%% 2D slices displayed 

% Reshape v image voxel data
v_vol = permute(reshape(v, Ny, Nx, Nz), [2 1 3]); % preserve order x,y,z

% Log compression (dB scale)
v_abs = abs(v_vol);
v_db = 20*log10(v_abs / max(v_abs(:)) + eps);   % eps avoids log(0)

% Display all Z-slices
figure;

nCols = 5; 
nRows = ceil(Nz / nCols);

for k = 1:Nz
    subplot(nRows, nCols, k);
    imagesc(x_range, y_range, v_db(:,:,k)');
    axis image off;
    caxis([-10 0]);   % 40 dB dynamic range
    title(sprintf('Z = %d', k));
end

colormap gray;
sgtitle(sprintf('LSQR Reconstruction (Iteration = %d)', maxIter));
colorbar;
%% Interactive Scroll through slices

figure;
for k = 1:Nz
    imagesc(x_range, y_range, abs(v_vol(:,:,k)'));
    axis image;
    colormap gray;
    colorbar;
    title(sprintf('Z-slice %d / %d', k, Nz));
    xlabel('x (mm)');
    ylabel('y (mm)');
    pause(0.3);  % adjust speed
end

%% log compressed bit
figure;
for k = 1:Nz
    % Extract and log-compress slice
    slice = abs(v_vol(:,:,k))';
    slice_db = 20*log10(slice / max(slice(:)) + eps);
    
    % Display
    imagesc(x_range*1e3, y_range*1e3, slice_db);
    axis image;
    colormap gray;
    caxis([-9 0]);          % 9 dB dynamic range
    colorbar;
    
    xlabel('x (mm)');
    ylabel('y (mm)');
    title(sprintf('Z-slice %d / %d', k, Nz));
    
    pause(0.3);               % adjust animation speed
end
%% 3D scatter plot
% voxel_positions = [x y z] in meters

% Normalised values
v_vol_norm = abs(v_vol);            % magnitude of complex signals
v_vol_norm = v_vol_norm / max(v_vol_norm(:));  % scale 0..1

% v_vol_norm(:) = normalized intensity values

figure;
scatter3(sensor_pos(:,1), ...   % x 
         sensor_pos(:,2), ...   % y 
         sensor_pos(:,3), ...   % z 
         100, ...                          % marker size
         v_vol_norm(:), ...               % color = intensity
         'filled');                       

xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
colormap(gray); colorbar;
title('Voxel-based 3D B-mode');
axis equal; view(3);

%% 3D scatter log compressed bit
% voxel_positions = [x y z] in meters

% Normalised values
v_vol_norm = abs(v_vol);            % magnitude of complex signals
v_vol_norm = v_vol_norm / max(v_vol_norm(:));  % scale 0..1

% Convert to dB scale (log-compression)
v_db = 20*log10(v_vol_norm(:) + eps);  % eps avoids log(0)

% Optional: clip dynamic range (e.g., 20 dB)
v_db_clip = max(v_db, max(v_db) - 10);

figure;
scatter3(voxel_positions(:,1), ...   % x 
         voxel_positions(:,2), ...   % y 
         voxel_positions(:,3), ...   % z 
         100, ...                       % marker size
         v_db_clip, ...                  % color = log-compressed intensity
         'filled');                     

xlabel('x (mm)'); ylabel('y (mm)'); zlabel('z (mm)');
colormap(gray); colorbar;
title('Voxel-based 3D B-mode (log-compressed dB)');
axis equal; view(3);
