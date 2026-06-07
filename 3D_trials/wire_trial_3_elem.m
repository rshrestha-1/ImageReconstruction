%%  Using kwave's data with 3 transducer elements:

load('sensor_data_filtered.ma');
load('voxel_map_x.ma');
load('voxel_map_y.ma');
load('voxel_map_z.ma');


param = [];
param.fc = 7.54e6;
param.bandwidth = 93;
param.width = 1.7000e-04;
param.height = 0.0050;
param.pitch = 2.0000e-04;
param.Nelements = 3;

param.xe = element_map_m(1,:);
param.ye = element_map_m(2,:);
param.elements = [param.xe; param.ye];
param.c = 1540;
param.fs = 30e6;

param.passive = true;


param.TXDELAY = zeros(1, param.Nelements);


IQb = das3(sensor_data_filtered,voxel_map_x, voxel_map_y, voxel_map_z, param);

I = 20 * log10(abs(IQb) / max(abs(IQb(:))));

figure
I(1:round(size(I,1)/2),1:round(size(I,2)/2),:) = NaN;
for k = [-40:10:-10 -5 -1]
      isosurface( ...
        voxel_map_x * 1e4, ...
        voxel_map_y * 1e4, ...
        voxel_map_z * 1e4, ...
        I, ...
        k ...
    );

end
colormap([1-hot;hot])
colorbar
box on, grid on
zlabel('[cm]')
title('PSF at (0,0,3) cm [dB]')
