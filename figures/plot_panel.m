function plot_panel(v_true, v_rec, scenario_name)

% Take central slice
Nz = size(v_true,3);
slice_idx = round(Nz/2);

gt_slice  = v_true(:,:,slice_idx);
rec_slice = v_rec(:,:,slice_idx);

% Difference (scaled for visibility)
diff_slice = abs(gt_slice - rec_slice);
diff_slice = 5 * diff_slice; % scaling factor

figure;

subplot(1,3,1)
imagesc(gt_slice)
axis image off
title([scenario_name ' - Ground Truth'])
colorbar

subplot(1,3,2)
imagesc(rec_slice)
axis image off
title('Reconstruction')
colorbar

subplot(1,3,3)
imagesc(diff_slice)
axis image off
title('Difference (x5)')
colorbar

colormap gray

end