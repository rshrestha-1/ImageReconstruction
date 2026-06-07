fs = 100e6;      % example: 100 MHz sampling frequency
f_low = 12e6;    % 12 MHz
f_high = 30e6;   % 30 MHz

% Design bandpass filter
Wn = [f_low f_high] / (fs/2);   % normalised frequencies
order = 4;                     % typical (4–6)

[b, a] = butter(order, Wn, 'bandpass');

% Filter per row (dimension 1)
u_filt = filtfilt(b, a, u);

% Test bandpass filter worked
freqz(b, a, 2048, fs)