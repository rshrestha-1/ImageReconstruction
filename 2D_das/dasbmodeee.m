%this code takes channel_data and finds BF_das through das function. then compares bmode images obtained from BF_das and das_data 


Ne = 128;
clear param

channel_data_2 = squeeze(channel_data(:,:,1,1));  % making it 2d

% better to have getparam()
param.fc = 7.5e6;
param.bandwidth = 93;
param.width = 170e-6;
param.height = 5e-3;
param.c = 1540;
param.fs = 40e6;
param.fnumber = 1.5;
param.Nelements = Ne;
param.pitch = 200e-6;

% Element positions
param.xe = ((0:Ne-1) - (Ne-1)/2) * param.pitch;
param.ze = zeros(1,Ne);
param.WasTransmitting = true(1,Ne);


param.TXDELAY = zeros(1,Ne);

x_axis = linspace(-0.02,0.02,248);
z_axis = linspace(0.01,0.05,248);
[x,z] = meshgrid(x_axis,z_axis);


[BF_das,param] = das(channel_data_2,x,z,param);
BF_das_1d = BF_das(:);

% Beamformed image
% subplot(122)
% pcolor(x*100,z*100,abs(BF_das).^.5)
% colormap(gray)
% xs = x;
% zs = z;
% hold on
% %plot(xs*100,zs*100,'rx')
% hold off
% title('Gamma-compressed image')
% xlabel('[cm]'), ylabel('[cm]')
% shading interp, axis equal ij tight

env_1 = abs(hilbert(BF_das));


bmode_1 = 20*log10(env_1 / max(env_1(:)));  
bmode_1 = bmode_1.'; %the image is more accurate with this
bmode_1(bmode_1 < -60) = -60;   

%  B-mode image
imagesc(x_axis*100, z_axis*100, bmode_1);
colormap(gray);
axis equal;
set(gca,'YDir','reverse');
title('B-mode image');
colorbar;

figure;
imagesc(abs(BF_das).'); %this shows where the scatterers are



%%trying with the das data

%pixel grid used in their das
x_axis_1 = linspace(-0.0132,0.0002,251);
z_axis_1 = linspace(0,0.0002,251);
[x_1,z_1] = meshgrid(x_axis_1,z_axis_1);


das_data_= squeeze(das_data(:,1,:));

env_2 = abs(das_data_);


bmode_2 = 20*log10(env_2 / max(env_2(:)));  
bmode_2 = bmode_2.';  %the image is more accurate with this


bmode_2(bmode_2 < -60) = -60;  

figure;
imagesc(x_axis*100, z_axis*100, bmode_2,[-50 0]);
colormap(gray);
shading interp;
axis equal ij tight;
set(gca,'YDir','reverse');

title('B-mode image');


figure;
imagesc(abs(das_data_).'); 

