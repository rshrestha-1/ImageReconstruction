close all; 
clear all; 
clc
load('H_exp_new.mat');
load('u_exp_new.mat'); %resin version
%% VOXEL MAP USED IN EXPERIMENTAL

Nx = 6; 
Ny = 6;
Nz = 6; 

N_voxels = Nx * Ny * Nz; 

x_range = linspace(-2.5e-3, 2.5e-3, Nx); 
y_range = linspace(-2.5e-3, 2.5e-3, Ny); 
z_range = linspace(14e-3, 19e-3, Nz); 

[x_grid, y_grid, z_grid] = ndgrid(x_range, y_range, z_range);


%% AUTOCONVOLUTION OF H


Ne = 58;

Nt = size(H_exp_new,1)/Ne;

Hf_3D = reshape(H_exp_new, [Nt, Ne, N_voxels]);

H_conv = zeros(Ne*Nt, N_voxels);

for v = 1:N_voxels
    for e = 1:Ne
      
        h_tx = Hf_3D(:, e, v);

        h_echo_full = conv(h_tx, h_tx);
        
        % truncate
        h_echo = h_echo_full(1:Nt);
        
        % Normalisation
        h_echo = h_echo / max(abs(h_echo) + eps);
        
        % put into initialised H
        row_start = (e-1)*Nt + 1;
        row_end   = e*Nt;
        H_conv(row_start:row_end, v) = h_echo; 
    end
end



%% Band-pass filtering 
% Fs = 30e6;                 
% u_filtered = zeros(size(u_exp_new));
% 
% bpFilt = designfilt('bandpassfir', 'FilterOrder', 100, 'CutoffFrequency1', 3e6,'CutoffFrequency2', 12e6,'SampleRate', Fs);
% 
% for k = 1:size(u_exp_new,2)
%     u_filtered(:,k) = filtfilt(bpFilt, u_exp_new(:,k)); 
% end

%% LSQR TO GET V


maxIter = 19;
tol = 1e-8;
lamdba = e-5;
lamdba_2 = e-10;

% v = twist_bpdn(H_conv, u_exp_new, lamdba,maxIter, tol);
%v = twist_elastic(H_conv, u_exp_new, lamdba,lamdba_2,maxIter, tol);
%v = fista_lasso(H_conv, u_exp_new, lambda, maxIter, tol);
%v = fista_elastic(H_conv, u_exp_new, lamdba,lamdba_2, maxIter, tol);

v = lsqr(H_conv, u_exp_new, tol, maxIter);

%% RESHAPE V
V = reshape(abs(v), Nx, Ny, Nz);

%% NEW FINE GRID

Nx_fine = 18;
Ny_fine = 18;
Nz_fine = 6;

x_fine = linspace(min(x_range), max(x_range), Nx_fine);
y_fine = linspace(min(y_range), max(y_range), Ny_fine);
z_fine = linspace(min(z_range), max(z_range), Nz_fine);

[xq, yq, zq] = ndgrid(x_fine, y_fine, z_fine);


%interpolate to the new grid
V = interpn(x_grid, y_grid, z_grid, V, xq, yq, zq, 'linear'); 

%% ANGULAR APPROACH

%spacing between x voxels
dx = 0.001;       
dy = 0.001;
      
lambda = 0.2e-3; 

%frequency index
k = 2*pi/lambda; 

% Fourier coordinates
fx = (-Nx_fine/2:Nx_fine/2-1)/(Nx_fine*dx);
fy = (-Ny_fine/2:Ny_fine/2-1)/(Ny_fine*dy);
[FX,FY] = meshgrid(fx, fy);

%Kz equation
kz = sqrt(k^2 - (2*pi*FX).^2 - (2*pi*FY).^2);

dz = 1e-3;  

%number of slices we want to propagate to
Nz_asp = 2;

V_asp = zeros(Nx_fine, Ny_fine, Nz_asp);
ind = 1;

%reference slice
slice0 = V(:,:,ind);

for zi = 1:Nz_asp

    F_slice = fftshift(fft2(slice0)); %take fourier transform

    F_slice_prop = F_slice .* exp(1i * kz * dz * zi); %find fourier transform of V for next slice

    V_asp(:,:,zi) = abs(ifft2(ifftshift(F_slice_prop))); %inverse fourier
end

%new z range
z_range_asp = z_range(ind) + (1:Nz_asp)*dz;

%% IMAGE AFTER ANGULAR APPRAOCH

V_env_1 = abs(V_asp);
V_env_1 = V_env_1 / max(V_env_1(:));

% Log compression (B-mode style)
dynamic_range_dB = 15;  
V_dB = 20*log10(V_env_1 + eps);  
% Clip dynamic range
V_dB(V_dB < -dynamic_range_dB) = -dynamic_range_dB;
% Plot
figure;
for k = 1:Nz_asp
    subplot(ceil(Nz_asp/4), 4, k); 
    imagesc(x_fine*1e3, y_fine*1e3, squeeze(V_dB(:,:,k))');
    axis image
    set(gca,'YDir','normal');   
    colormap(gray)
    clim([-dynamic_range_dB 0]);
    datacursormode on
    title(sprintf('z = %.2f mm, at slice = %d', z_range_asp(k)*1e3, k));
    
    if k == 1
        xlabel('x (mm)')
        ylabel('y (mm)')
    end
end

sgtitle(['Reconstruction with ' num2str(maxIter) ' iterations and ' num2str(dynamic_range_dB) 'dB' ])
%% ORIGINAL IMAGE BEFORE ANGULAR


V_env_1 = abs(V);
V_env_1 = V_env_1 / max(V_env_1(:));

% Log compression (B-mode style)
dynamic_range_dB = 10;  
V_dB = 20*log10(V_env_1 + eps);  
% Clip dynamic range
V_dB(V_dB < -dynamic_range_dB) = -dynamic_range_dB;
% Plot
figure;
for k = 1:Nz_fine
    subplot(ceil(Nz/4), 4, k); 
    imagesc(x_fine*1e3, y_fine*1e3, squeeze(V_dB(:,:,k))');
    axis image
    set(gca,'YDir','normal');   
    colormap(gray)
    clim([-dynamic_range_dB 0]);
    datacursormode on
    title(sprintf('z = %.2f mm, at slice = %d', z_fine(k)*1e3, k));
    
    if k == 1
        xlabel('x (mm)')
        ylabel('y (mm)')
    end
end

sgtitle(['Reconstruction with ' num2str(maxIter) ' iterations and ' num2str(dynamic_range_dB) 'dB' ])

%% Check if H has imaginary values

if ~isreal(H_exp_new)
    disp('H contains imaginary data')
    
    % Count elements with imaginary part
    numComplex = sum(abs(imag(H_exp_new(:))) > 1e-12); % use threshold to catch tiny floating-point values
    
    % Correct way to print formatted text
    fprintf('Number of elements with imaginary part: %d\n', numComplex)
else
    disp('H is purely real')
end
