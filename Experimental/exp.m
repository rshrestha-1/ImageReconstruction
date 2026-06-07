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
%%
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
