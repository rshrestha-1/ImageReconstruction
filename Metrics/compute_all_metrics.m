function results = compute_all_metrics(v_true, v_rec, labels)

% Ensure vectors
v_true = double(v_true);
v_rec  = double(v_rec);

%% MSE 
results.MSE = mean((v_true(:) - v_rec(:)).^2);

%% PSNR
max_val = max(v_true(:));
results.PSNR = 10 * log10(max_val^2 / results.MSE);

%% SSIM (slice-wise)
Nz = size(v_true,3);
ssim_vals = zeros(Nz,1);

for k = 1:Nz
    ssim_vals(k) = ssim(v_rec(:,:,k), v_true(:,:,k));
end

results.SSIM = mean(ssim_vals);

%% Region masks
blood_mask = labels == "blood";
lipid_mask = labels == "lipid";
calc_mask  = labels == "calcified";

%% SNR (plaque vs blood noise) 
plaque_mask = lipid_mask | calc_mask;

mu_signal = mean(v_rec(plaque_mask));
sigma_noise = std(v_rec(blood_mask));

results.SNR = mu_signal / sigma_noise;

%% CNR (plaque vs blood) 
mu_blood  = mean(v_rec(blood_mask));
mu_plaque = mean(v_rec(plaque_mask));
sigma_all = std(v_rec(:));

results.CNR = abs(mu_plaque - mu_blood) / sigma_all;

%% Lipid contrast 
if any(lipid_mask(:))
    mu_lipid = mean(v_rec(lipid_mask));
    results.LipidContrast = abs(mu_lipid - mu_blood) / ...
        (mu_lipid + mu_blood);
else
    results.LipidContrast = NaN;
end

%% Calcified contrast
if any(calc_mask(:))
    mu_calc = mean(v_rec(calc_mask));
    results.CalcifiedContrast = abs(mu_calc - mu_blood) / ...
        (mu_calc + mu_blood);
else
    results.CalcifiedContrast = NaN;
end

end
