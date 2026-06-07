% Adjoint operator (transpose) 
H_adj = H';   % H*  (adjoint)
HTH = H_adj * H; % HTH
HTH_inv = inv(HTH); % HTH_inv
H_pinv = HTH_inv * H_adj; % Pseudo-inverse H

% Equivalent Built-in MATLAB Function
H_builtin_pinv = pinv(H)

% Use for u to calculate v
v = H_pinv * u;
v_builtin_pinv = H_builtin_pinv * u
