id="setup"
fs = 30e6;          % sampling frequency
f_low = 5e6;
f_high = 10e6;

lambda1 = 0.01;     % L1 (sparsity)
lambda2 = 0.001;    % L2 (stability)
alpha = 1e-3;       % step size
nIter = 100;
%% 
plot(cost)
xlabel('Iteration')
ylabel('Cost')
title('Convergence of FISTA')


%% smooth u
u_filt = smoothdata(u, 'gaussian', 5);
%% initiakise lsqr light regularisation
x = lsqr(Hn, b, 1e-6, 10);  % low iterations only

%% fft u

fs = 30e6; % 30 MHz
N = length(u);

U = fft(u);
f = (0:length(u)-1)*(fs/length(u));

figure;
plot(f/1e6, abs(U))
xlabel('Frequency (MHz)')
ylabel('Magnitude')
title('FFT of u')
xlim([0 fs/2]/1e6) % only show positive frequencies

% FFT
U = fft(u);
P = abs(U).^2; % Power spectrum

% Only use positive frequencies
half = 1:floor(N/2);
f_half = f(half);
P_half = P(half);

% Normalise power
P_norm = P_half / max(P_half);

% Index where power > threshold (e.g. -20 dB ≈ 0.01)
threshold = 0.01;

idx = find(P_norm > threshold);

f_low = f_half(min(idx));
f_high = f_half(max(idx));

fprintf('Optimal band: %.2f MHz to %.2f MHz\n', f_low/1e6, f_high/1e6);
%
figure;
plot(f_half/1e6, 10*log10(P_half/max(P_half)));
hold on;
xline(f_low/1e6, 'r--');
xline(f_high/1e6, 'r--');

xlabel('Frequency (MHz)');
ylabel('Power (dB)');
title('Detected Bandwidth');
grid on;
%% tik ridge
figure;


% L2 sweep 1e-4, 1e-3, 1e-2
lambda2 = 1e-3;

% Tik ridge with lsqr
A = [Hn; sqrt(lambda2)*eye(size(Hn,2))];
b_aug = [u; zeros(size(Hn,2),1)];

v_l2 = lsqr(A, b_aug, 1e-6, 50);
plot(abs(v_l2));

figure;

% L1 soft thresholding test 0.01, 0.03, 0.05, 0.1
lambda1 = 0.05 * max(abs(v_l2));

v_sparse = sign(v_l2) .* max(abs(v_l2) - lambda1, 0);
plot(abs(v_sparse));
%% 
load('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/sensor_data_test.mat')
u = sensor_data_ds;
%% 
%% reshape H
H_reshaped = reshape(H, 10000, 90, 216);
H_cut = H_reshaped(1:8192, 1:58, :);

H_exp_new = reshape(H_cut, 8192*58, 216);
