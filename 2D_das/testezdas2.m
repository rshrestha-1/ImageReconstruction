Ne = 128;
clear param

channel_data_2 = squeeze(channel_data(:,:,1,1));  

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

param.xe = ((0:Ne-1) - (Ne-1)/2) * param.pitch;
param.ze = zeros(1,Ne);
param.WasTransmitting = true(1,Ne);

% temp delays
param.TXDELAY = zeros(1,Ne);

x_axis = linspace(-0.02,0.02,248);
z_axis = linspace(0.01,0.05,248);
[x,z] = meshgrid(x_axis,z_axis);


[BF_ezdas,param] = ezdas2(channel_data_2,x,z,[0,0],param);
BF_ezdas_1d = BF_ezdas(:);



env_1 = abs(hilbert(BF_ezdas));


bmode_1 = 20*log10(env_1 / max(env_1(:)));  
bmode_1 = bmode_1.'; 
bmode_1(bmode_1 < -60) = -60;   

%  B-mode image
imagesc(x_axis*100, z_axis*100, bmode_1);
colormap(gray);
axis equal;
set(gca,'YDir','reverse');
title('B-mode image');
colorbar;

figure;
imagesc(abs(BF_ezdas).'); %this shows where the scatterers are

