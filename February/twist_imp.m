Nx = 38;
Ny = 17;
Nz = 8;
N_voxels = Nx*Ny*Nz;

true_idx = 513;
all_idx  = 1:N_voxels;
bg_idx   = setdiff(all_idx, true_idx);

lambdas = logspace(-5, 0, 20);
contrast = zeros(size(lambdas));

L = norm(H_f)^2;  % squared spectral norm used in twist
tau = 1 / L;

A  = @(x) H_f * x; %
AT = @(x) H_f' * x;


for i = 1:length(lambdas)


    
    v = twist_bpdn(A, AT, u, lambdas(i), 100,tau);

    v = v(1:N_voxels);  % ensure correct length

    I_signal = abs(v(true_idx));
    I_bg     = mean(abs(v(bg_idx)));

    contrast(i) = 20*log10(I_signal / I_bg);

end

[best_contrast, idx] = max(contrast);
best_lambda = lambdas(idx);

%% 
v_best = twist_bpdn(A, AT, u, best_lambda, 100,tau);
