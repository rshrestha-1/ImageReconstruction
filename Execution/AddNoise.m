%% Add 10% signal peak as white Gaussian noise

peak_signal = max(abs(H(:)));
sigma = 0.01*peak_signal;
noise = sigma*randn(size(H));
H = H+noise;
