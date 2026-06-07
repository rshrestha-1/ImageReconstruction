%% Bmode slices (2D)

%% Solution via LSQR

maxIter = 18;
tol = 1e-6;

v = lsqr(H,u,tol,maxIter);

%% Filter image v

fs = 30e6;         
f_l  = 3e6;      
f_h = 12e6;      

bpFilt = designfilt('bandpassiir', 'FilterOrder', 6, 'HalfPowerFrequency1', f_l, 'HalfPowerFrequency2', f_h, 'SampleRate', fs);

v_filt = filtfilt(bpFilt, v); 

%% Display image v as bmode

% B-mode image
env = abs(hilbert(v_filt)); % amplitude envelope
env_norm = env ./ max(env(:)); % normalised envelope

dynamicRange = 9; % dB (typical: 40–80 dB)

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