function err = localisation_error(labels_true, labels_rec, type)

% centroid of true
[idx_x, idx_y, idx_z] = ind2sub(size(labels_true), find(labels_true == type));
centroid_true = [mean(idx_x), mean(idx_y), mean(idx_z)];

% centroid of reconstructed
[idx_x, idx_y, idx_z] = ind2sub(size(labels_rec), find(labels_rec == type));
centroid_rec = [mean(idx_x), mean(idx_y), mean(idx_z)];

err = norm(centroid_true - centroid_rec);

end