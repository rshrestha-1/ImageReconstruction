function [SNR_val, CNR_val] = compute_snr_cnr(v_rec, labels)

% Example: blood vs plaque
blood = v_rec(labels == "blood");
plaque = v_rec(labels ~= "blood");

mu_signal = mean(plaque);
sigma_noise = std(blood);

% SNR
SNR_val = mu_signal / sigma_noise;

% CNR
mu1 = mean(plaque);
mu2 = mean(blood);
sigma = std(v_rec(:));

CNR_val = abs(mu1 - mu2) / sigma;

end