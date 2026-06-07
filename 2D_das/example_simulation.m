Ne = 128;
clear param

% Remove angles 5–14 index
new_data = channel_data;
new_data(:,:,:,5:14) = [];

f_low  = 2e6;
f_high = 6e6;             
param.fs = settings.acquisition.fs;

new_data = double(new_data);

bpFilt = designfilt('bandpassiir', ...
    'FilterOrder', 6, ...
    'HalfPowerFrequency1', f_low, ...
    'HalfPowerFrequency2', f_high, ...
    'SampleRate', param.fs);

filtered_data = zeros(size(new_data));

for ia = 1:size(new_data,4)
    for ch = 1:size(new_data,2)
        sig = new_data(:,ch,1,ia);
        filtered_data(:,ch,1,ia) = filtfilt(bpFilt, sig);
    end
end

new_data = filtered_data;
clear filtered_data        


%new_data = sensor_data_2;
angles_kept = settings.acquisition.angles;
angles_kept(5:14) = [];

% apply Lens correction
lens_corr_mm = 1.183;
lens_corr_m  = lens_corr_mm * 1e-3;

param.c = 1540;
% param.t0 = lens_corr_m / param.c;
param.t0 = settings.acquisition.t0(1);
% Parameters 
param.fc = 7.5e6;
param.bandwidth = 93;
param.width = 170e-6;
param.height = 5e-3;
param.fs = 15e6;
param.fnumber = 1.5;
param.Nelements = Ne;
param.pitch = 200e-6;

param.xe = ((0:Ne-1) - (Ne-1)/2) * param.pitch;
param.ze = zeros(1,Ne);
param.WasTransmitting = true(1,Ne);


x_axis_1 = linspace(-13e-3, 13e-3, 251);
z_axis_1 = linspace(0, 50e-3, 251);  
[x_1, z_1] = meshgrid(x_axis_1, z_axis_1);


%  Beamform 
Nangles = size(new_data,4);
BF_compounded = zeros(size(x_1));

for ia = 1:Nangles

    SIG = squeeze(new_data(:,:,1,ia));
    theta = angles_kept(ia);

    tx_delay = (param.xe * sin(theta)) / param.c;
    tx_delay = tx_delay - min(tx_delay);

    param.TXdelay = tx_delay;   

    [BF, param] = das(SIG, x_1, z_1, param);

    BF_compounded = BF_compounded + BF;
end

env_2 = abs(BF_compounded);


bmode_2 = 20*log10(env_2 / max(env_2(:)));  


bmode_2(bmode_2 < -80) = -80;  


figure;
%imagesc(x_axis_1, z_axis_1, bmode_2,[-80 0]);

imagesc(x_axis_1*100, z_axis_1*100, bmode_2,[-50 0]);
colormap(gray(256));
shading interp;
axis equal ij tight;
set(gca,'YDir','reverse');

title('B-mode image');
