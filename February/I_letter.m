close all; 
clear all; 
clc

load('udata_logo.mat');
load('H_f.mat')
u = udata_logo;
n_elements = 3;
time_samples = 1830;

%% 
H_conv = ifft2( fft2(H_f) .* fft2(H_f) );
H_conv = real(H_conv); 
%% 
tol = 1e-6;
maxit = 18;
%%
%% 

v = lsqr(H_conv, u, tol, maxit);
%% 

v = H_conv'*u;
%% 
Nx = 16; 
Ny = 16; 
Nz = 20; 
%% 

V = reshape(v(1:Nx*Ny*Nz), Ny, Nx, Nz);
x_range = linspace(-4.4e-3, 4.4e-3, Nx); 
y_range = linspace(-4.4e-3, 4.4e-3, Ny); 
z_range = linspace(1e-3, 4e-3, size(V,3)); 
[xi, yi, zi] = meshgrid(x_range, y_range, z_range);

%% 
figure;
slice(xi, yi, zi, V, [], [], z_range);
shading interp;
colorbar;
view(3);
%% 
figure;

for k = 1:Nz
    subplot(1, Nz, k)
    imagesc(V(:,:,k))
    axis image off
end

colormap 
