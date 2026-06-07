%% Autoconvolved H fully
Ne = 3;        % number of elements
Nt = 1830;     % number of time samples per element
Nv = 5168;     % number of voxels

Hf_3D = reshape(H_f, [Nt, Ne, Nv]);

H = zeros(Ne*Nt, Nv);

for v = 1:Nv
    
    for e = 1:Ne
      
        h_tx = Hf_3D(:, e, v);
        
        % autoconvolution
        h_echo_full = conv(h_tx, h_tx);
        
        % truncate
        h_echo = h_echo_full(1:Nt);
        
        % Normalisation
        % h_echo = h_echo / max(abs(h_echo) + eps);
        
        % put into initialised H
        row_start = (e-1)*Nt + 1;
        row_end   = e*Nt;
        
        H(row_start:row_end, v) = h_echo;
        
    end
    
end

