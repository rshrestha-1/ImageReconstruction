%% Optimising lambda

% Define ROI for letter
ROI_letter = false(19,19,22); % 3D array of volume Nx Ny Nz = filled with false
ROI_letter(9:10,9:10,1:2) = true; % region containing letters (becomes true)

% Define ROI for background
ROI_background = false(19,19,22); % same imaging volume set to false
ROI_background(18:19,14:19,16:19) = true;% set corner of volume with no object for background intensity

% Return index of positions where true
roi_letter_idx = find(ROI_letter);
roi_bg_idx = find(ROI_background);

% Compute average letter intensity
I_letter = mean(abs(v(roi_letter_idx)));
I_bg = mean(abs(v(roi_bg_idx)));

% Sweep lambda values
lambdas = logspace(-5,0,20);

for i = 1:length(lambdas)
    % Apply TwIST at the lambda i for 100 iterations
    v = twist_bpdn(H,H',u,lambdas(i),100,1);
    
    % Compute average letter intensity
    I_letter = mean(abs(v(roi_letter_idx)));
    I_bg = mean(abs(v(roi_bg_idx)));

    % Compute contrast and store per lambda
    contrast(i) = 20*log10(I_letter/I_bg);
end

% Identify best contrast ratio and lambda value to obtain it
[best_contrast,idx] = max(contrast);
best_lambda = lambdas(idx);
