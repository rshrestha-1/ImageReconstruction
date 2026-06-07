function [bfSIG,M] = ezdas(SIG,x,y,vsource,param)

% EZDAS Easier DAS
% 
% BFSIG = EZDAS(SIG,X,Z,VSOURCE,PARAM) beamforms 
% the RF or I/Q signals stored in the array SIG,
% and returns the beamformed signals BFSIG. The
% signals are beamformed at the points specifed by X and X.
% 
% 1) SIG must be a 2D array. The first dimension
% (i.e. each column) corresponds to single RF
% or I/Q signals over (fast-) time, with the
% 1st column corresponding to the 1st element.
% 2 VSOURCE contains the coordinates [x0,z0] of
% the virtual point source. Use large x0,z0 for the
% plane waves.
% 3) PARAM is a structure that contains the 
% paramter values required for the delaty-and-sum
% 
% Note: SIG must be complex for I/Q data
% (i.e. SIG = complex(I,Q) = I + li*Q).
% 
% [~,M] = EZDAS(...) also returns the DAS matrix)
% EZDAS uses a linear interpolation to generate the DAS matrix

% PARAM Fields:
% PARAM.fs
% PARAM.pitch
% PARAM.fc
% PARAM.c
% PARAM.fnumber

% INPUTS:
% SIG: nSamples × nChannels matrix
% Each column is the RF/IQ signal from 1 array element
%
% x, z: coordinates of the imaging points (grids). Can be 2-D matrices.
% X and Z must match the desired beamforming image size.
%
% vsource: [x0, z0] coordinates of virtual source (plane wave -> very far)

% OUTPUTS:
% bfSIG: beamformed image (same size as x and z)
% M: sparse DAS matrix that performs the operation bfSIG = M * SIG(:)

% ADD DELAY FROM MASK

% Reshape coordinates
siz0 = size(x); % size of the output image grid
[n1,nc] = size(SIG);
% convert grid to column vectors
x = x(:); z = z(:);

% ULA (uniform linear array):
% x coordinates of the elements
xe = ((0:nc-1)-(nc-1)/2)*param.pitch;
L = xe(end)-xe(1); % length of the array

% coordinates  of the virtual source
x0 = vsource(1); z0 = vsource(2); y0 = 0

% transmit & receive distances

% 1) TRANSMIT DISTANCE: voxel <- virtual source
% dTX = sqrt( (x - x0)^2 + (z - z0)^2 + (y - y0)^2 )

dTX = hypot(x - x0, z - z0, y - y0);
hypot((abs(x0)-L/2)*(abs(x0)>L/2),z0,y0);

% 2) RECEIVE DISTANCE: voxel -> each array element
% dRX(m) = sqrt( (x - xe(m))^2 + z^2 + y^2)

dRX = hypot(x-xe,z,y);

% 3) TOTAL TWO-WAY TRAVEL TIME
% tau = (dTX + dRX) / c

tau = (dTX+dRX)/param.c; % two-way travel times

% fast-time indices
idxt = tau*param.fs + 1; % Continuous sample index (1-based)

% boolean vectors
I = idxt>=1 & idxt<=n1-1; % Must be inside valid time range
% Receive aperture condition:
Iaperture = abs(x-xe)<=(z/2/param.fnumber); % |x - xe| <= z / (2*fnumber)
I = I&Iaperture; % Final boolean mask (only valid + inside aperture)

% linear indices
idx = idxt + (0:nc-1)*n1; % Flatten everything that passes mask
idx = idx(I); % continuous index
idxf = floor(idx); % integer part
idx = idxf-idx; % interpolation fraction

% DAS matrix
[i,~] = find(I);
s = [idx+1;-idx]; % (for linear interpolation)
if ~isreal(SIG) % if IQ: phase rotations
    % Phase shift for IQ: exp(j 2π fc τ)
    s = s.*exp(2i*pi*param.fc*[tau(I);tau(I)]);
end

M = sparse([i;i],[idxf;idxf+1],s,numel(x),n1*nc);

% DAS beamforming
bfSIG = reshape(M*SIG(:),siz0);
