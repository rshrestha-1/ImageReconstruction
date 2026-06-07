function metrics = compute_metrics(v_true, v_rec)

% MSE
metrics.MSE = mean((v_true(:) - v_rec(:)).^2);

% PSNR
max_val = max(v_true(:));
metrics.PSNR = 10 * log10(max_val^2 / metrics.MSE);

% SSIM (slice-wise average)
Nz = size(v_true,3);
ssim_vals = zeros(Nz,1);

for k = 1:Nz
    ssim_vals(k) = ssim(v_rec(:,:,k), v_true(:,:,k));
end

metrics.SSIM = mean(ssim_vals);

end