%% Visualise voxel sensitity
sens = reshape(sum(abs(H),1),Nx,Ny,Nz);
figure;
montage(sens,'DisplayRange',[])
%% Column Normalisation
col_norms = sqrt(sum(H.^2,1));
Hn = H ./ col_norms;

sensn = reshape(sum(abs(Hn),1),Nx,Ny,Nz);

figure;
montage(sensn,'DisplayRange',[])
title('Voxel sensitivity after column normalisation')
%% Test column correlation

rand_vox = randperm(Nx*Ny*Nz,400);

C = corrcoef(H(:,rand_vox));
figure;
imagesc(C)
colorbar

C_no_diag = C - diag(diag(C));

max_corr = max(abs(C_no_diag(:)));
mean_corr = mean(abs(C_no_diag(:)));

fprintf('Max corr: %.3f\n', max_corr);
fprintf('Mean corr: %.3f\n', mean_corr);
