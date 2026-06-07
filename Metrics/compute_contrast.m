function contrast = compute_contrast(v_rec, labels, type1, type2)

region1 = v_rec(labels == type1);
region2 = v_rec(labels == type2);

contrast = abs(mean(region1) - mean(region2)) / ...
           (mean(region1) + mean(region2));

end