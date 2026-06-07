%this code works with 128 elements, and compares the data we would get from das function with the das data by haviinng correspondong b mode images

Ne = 128;
clear param


channel_data_2 = squeeze(channel_data(:,:,1,1));  % our fake data

% better to have getparam()
param.fc = 7.5e6;
param.bandwidth = 93;
param.width = 170e-6;
param.height = 5e-3;
param.c = 1540;
param.fs = 15000000;
param.fnumber = 1.5;
param.Nelements = Ne;
param.pitch = 200e-6;

% Element positions
param.xe = ((0:Ne-1) - (Ne-1)/2) * param.pitch;
param.ze = zeros(1,Ne);
param.WasTransmitting = true(1,Ne);


zf = 0.03; % this is just to test
xe = ((0:Ne-1)-(Ne-1)/2)*param.pitch;
param.TXDELAY = sqrt(param.xe.^2 + zf^2)/param.c; 
param.TXDELAY = param.TXDELAY - min(param.TXDELAY);




x_axis = linspace(-0.02,0.02,243);
z_axis = linspace(0.01,0.05,248);
[x,z] = meshgrid(x_axis,z_axis);


[BF_das,param] = das(channel_data_2,x,z,param);
BF_das_1d = BF_das(:);


%bmodeee

env_1 = abs(hilbert(BF_das));


bmode_1 = 20*log10(env_1 / max(env_1(:)));  
%bmode_1 = bmode_1.'; 
bmode_1(bmode_1 < -60) = -60;   


imagesc(x_axis*100, z_axis*100, bmode_1);
colormap(gray);
axis equal;
set(gca,'YDir','reverse');
title('B-mode image');
colorbar;


 % to see scatterers clearly
figure;
imagesc(x_axis*100, z_axis*100, abs(BF_das));  
colormap(gray);
axis equal;
set(gca,'YDir','reverse');  % depth increases downward

title('scatterers for data from data function');
colorbar;


%trying with the das data

%%bmode
x_axis_1 = linspace(-0.0132,0.0002,251);
z_axis_1 = linspace(0,0.0002,251);
[x_1,z_1] = meshgrid(x_axis_1,z_axis_1);


das_data_= squeeze(das_data(:,1,:));
das_data_rev = das_data_.';  % transpose


env_2 = abs(das_data_);


bmode_2 = 20*log10(env_2 / max(env_2(:)));  
bmode_2 = bmode_2.';  


bmode_2(bmode_2 < -60) = -60;  


figure;
imagesc(x_axis*100, z_axis*100, bmode_2,[-50 0]); %try x_axis_1 and z_axis_1 so it can work
colormap(gray);
shading interp;
axis equal ij tight;
set(gca,'YDir','reverse');


%scatterers again to compare

figure;
imagesc(x_axis*100, z_axis*100, abs(das_data_rev));  % do NOT transpose
colormap(gray);
axis equal;
set(gca,'YDir','reverse');  % depth increases downward

title('Scatteres for the das data');
colorbar;
