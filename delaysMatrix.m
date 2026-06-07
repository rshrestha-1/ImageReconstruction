
Ne = 192;             
pitch  = 200e-6;        

% coordinates of elements
x_i = zeros(1, Ne);
y_i = zeros(1, Ne);        
z_i = zeros(1,Ne);


%calculating coordinates using the "so you think you can das"
for idx = 1:Ne
    x_i(idx) = (pitch/2) * (2*idx - Ne - 1);
    y_i(idx) = 0;
    z_i(idx) = 0;
end

param = getparam('L12-3V');

param.elements = [x_i; y_i;z_i];


% make an image volume

Nx = 10; Ny = 10; Nz = 10;
numVoxels = Nx * Ny * Nz;

ijk = zeros(numVoxels, 3);

x_range = linspace(-5e-3, 5e-3, Nx);
y_range = linspace(0, 10e-3, Ny);
z_range = linspace(0, 20e-3, Nz);

[x, y, z] = ndgrid(x_range, y_range, z_range);


for n = 1:numVoxels
    [i, j, k] = ind2sub([Nx Ny Nz], n);
    ijk(n, :) = [i, j, k];
end

nSamples  = 1 % as a starting point


%delay matrix
delays_matrix = zeros(numVoxels,Ne);


for pointIdx = 1:numel(x)
    a = ijk(pointIdx,1);
    b = ijk(pointIdx,2);
    c = ijk(pointIdx,3);

    delays_matrix(pointIdx, :) = txdelay3(x(a,b,c), y(a,b,c), z(a,b,c), param);  %creates sparse delay matrix: [1000,192]

end  

param.fs = 9;

t = (0:nSamples-1)'/param.fs;

pulse = sin(2*pi*param.fc*t) .* exp(-(t-5e-6).^2/(0.3e-6)^2);

u = repmat(pulse,1,nChannels) + 0.05*randn(nSamples,Ne);

u = u(:); %this makes u a one column arrray [192x1]

%linear model, v is a [1000x1] matrix

v = delays_matrix*u;
