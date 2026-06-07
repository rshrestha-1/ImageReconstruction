

%L12-3V parameters
param = getparam('L12-3V');
% random RF signal
param.fs = 9;
param.c = 999;
param.fnumber = 0.01;
nSamples  = 1024;    % very random number: 100 samples per element
nChannels = 64;      % number of elements (typical L12 probe)

t = (0:nSamples-1)'/param.fs;

pulse = sin(2*pi*param.fc*t) .* exp(-(t-5e-6).^2/(0.3e-6)^2);

u = repmat(pulse,1,nChannels) + 0.05*randn(nSamples,nChannels);


% temporary image volume

x_vec = linspace(-10e-3,10e-3,128);
y_vec = linspace(-10e-3,10e-3,128);
z_vec = linspace(5e-3,40e-3,128);

[x, y, z] = meshgrid(x_vec, y_vec, z_vec);


[beamformed_signals,Das_matrix] = ezdas(u,x,y,z,[0,0],param);





