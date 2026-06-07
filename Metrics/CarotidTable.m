scenario_names = ["Healthy","Lipid","Fibrous","Calcified","Mixed"];

num_cases = length(scenario_names);

MSE = zeros(num_cases,1);
SSIM = zeros(num_cases,1);
SNR  = zeros(num_cases,1);
CNR  = zeros(num_cases,1);
LipC = zeros(num_cases,1);
CalC = zeros(num_cases,1);

for i = 1:num_cases

    % Generate phantom per scenario
    [v_true_3D, labels] = carotid_phantom(Nx,Ny,Nz, scenario_names(i));
    v_true = v_true_3D(:);

    % Forward model
    u = Hn * v_true;

    % % Add noise
    % SNR_dB = 20;
    % signal_power = mean(u.^2);
    % noise_power = signal_power / (10^(SNR_dB/10));
    % u_noisy = u + sqrt(noise_power)*randn(size(u));

    % Reconstruction
    lambda = 1e-3;
    v_rec = (Hn'*Hn + lambda*eye(size(Hn,2))) \ (Hn'*u_noisy);
    v_rec_3D = reshape(v_rec, size(v_true_3D));

    % Metrics
    res = compute_all_metrics(v_true_3D, v_rec_3D, labels);

    MSE(i)  = res.MSE;
    SSIM(i) = res.SSIM;
    SNR(i)  = res.SNR;
    CNR(i)  = res.CNR;
    LipC(i) = res.LipidContrast;
    CalC(i) = res.CalcifiedContrast;

end

% Table
ResultsTable = table(scenario_names', MSE, SSIM, SNR, CNR, LipC, CalC, ...
    'VariableNames', {'Scenario','MSE','SSIM','SNR','CNR','LipidContrast','CalcifiedContrast'});

disp(ResultsTable)