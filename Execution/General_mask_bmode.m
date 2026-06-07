%% Matlab Code
close all; clear all; clc
 
load('sensor_data.mat'); % Load the data

%% Load in u and H matrices

u = sensor_data_2;

H = eye(1525); % 1525 × 1525 (Identity Matrix as placeholder)

%% Via inverse

H_inv = inv(H); % remains Identity currently

v = H_inv * u; % size = 1525 × 3

%% Via adjoint

% Adjoint estimation v_est = H^H * u
v = H' * u;  % adjoint = transpose for real matrices
%filtering


fs = 30e6;         
f_l  = 3e6;      
f_h = 12e6;      

bpFilt = designfilt('bandpassiir', 'FilterOrder', 6, 'HalfPowerFrequency1', f_l, 'HalfPowerFrequency2', f_h, 'SampleRate', fs);

v_filt = filtfilt(bpFilt, v); 
%% Display image v as bmode

% B-mode image
env = abs(hilbert(v_filt)); % amplitude envelope
env_norm = env ./ max(env(:)); % normalised envelope

dynamicRange = 40; % dB (typical: 40–80 dB)

bmode = 20*log10(env_norm + eps); % log compression
bmode = bmode + dynamicRange;
bmode(bmode < 0) = 0;

figure
imagesc(bmode)
colormap(gray)
colorbar
clim([0 dynamicRange])
xlabel('Scan line')
ylabel('Depth')
title('B-mode Image')
