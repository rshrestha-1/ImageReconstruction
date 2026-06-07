
%% -------------------------------------------------------------------------
% Example:
% --------
% Plane-wave ultrasound simulation using a 3D matrix array
% and reconstruction of the 3D point-spread function (PSF)
% at 3 cm depth using delay-and-sum beamforming
% -------------------------------------------------------------------------

% Initialise parameter structure (used by SIMUS3, RF2IQ, DAS3)
param = [];

% -------------------------------------------------------------------------
% TRANSDUCER PARAMETERS
% -------------------------------------------------------------------------

% Central frequency of the transducer (3 MHz)
param.fc = 3e6;

% Fractional bandwidth of the transducer in percent
param.bandwidth = 70;

% Width of each transducer element [m]
param.width = 250e-6;

% Height of each transducer element [m]
param.height = 250e-6;

% -------------------------------------------------------------------------
% MATRIX ARRAY GEOMETRY
% -------------------------------------------------------------------------

% Pitch (centre-to-centre spacing) between elements [m]
pitch = 300e-6;

% Create a 32x32 grid of element positions in x and y
% The grid is centred around (0,0)
[xe, ye] = meshgrid(((1:32) - 16.5) * pitch);

% Store element positions as a 2 x Nelements matrix:
%   row 1 -> x positions
%   row 2 -> y positions
param.elements = [xe(:).'; ye(:).'];

% -------------------------------------------------------------------------
% POINT SCATTERER (PSF) LOCATION
% -------------------------------------------------------------------------

% Position of the point scatterer [m]
x0 = 0;        % lateral x
y0 = 0;        % lateral y
z0 = 3e-2;     % depth (3 cm)

% -------------------------------------------------------------------------
% TRANSMIT SETTINGS
% -------------------------------------------------------------------------

% Transmit time delays for plane-wave transmission
% (all zeros = plane wave)
dels = zeros(1, 1024);   % 32 x 32 = 1024 elements

% -------------------------------------------------------------------------
% RF DATA SIMULATION
% -------------------------------------------------------------------------

% Simulate raw RF channel data received by each element
% using the SIMUS3 ultrasound simulation
% RF: raw RF signals
% param: updated parameter structure
[RF, param] = simus3(x0, y0, z0, 1, dels, param);

% -------------------------------------------------------------------------
% RF → IQ DEMODULATION
% -------------------------------------------------------------------------

% Convert RF signals to complex baseband IQ data
% This performs:
%   - bandpass filtering
%   - quadrature demodulation
%   - envelope-ready representation
IQ = rf2iq(RF, param);

% -------------------------------------------------------------------------
% 3D VOXEL GRID FOR BEAMFORMING
% -------------------------------------------------------------------------

% Acoustic wavelength [m]
lambda = 1540 / param.fc;   % c / f

% Define 3D reconstruction grid:
%   x: lateral dimension
%   y: elevation dimension
%   z: depth dimension
[xi, yi, zi] = meshgrid( ...
    -2e-2 : lambda : 2e-2, ...           % x range (±2 cm)
    -2e-2 : lambda : 2e-2, ...           % y range (±2 cm)
     2.5e-2 : lambda/2 : 3.1e-2 );       % z range around 3 cm

% -------------------------------------------------------------------------
% 3D DELAY-AND-SUM BEAMFORMING
% -------------------------------------------------------------------------

% Perform 3D DAS beamforming:
%   - applies geometric delays
%   - sums contributions from all elements
% Result is a complex-valued 3D image
IQb = das3(IQ, xi, yi, zi, dels, param);

% -------------------------------------------------------------------------
% B-MODE FORMATION (ENVELOPE + LOG COMPRESSION)
% -------------------------------------------------------------------------

% Compute the B-mode image:
%   - abs(IQb): envelope detection
%   - normalisation by maximum value
%   - 20*log10: logarithmic compression (dB scale)
I = 20 * log10(abs(IQb) / max(abs(IQb(:))));

% -------------------------------------------------------------------------
% PSF VISUALISATION USING ISOSURFACES
% -------------------------------------------------------------------------

% Create a new figure
figure

% Mask one quadrant of the volume to avoid symmetric duplicates
I(1:round(size(I,1)/2), ...
  1:round(size(I,2)/2), :) = NaN;

% Plot isosurfaces at different dB levels
% These surfaces show the 3D PSF shape and sidelobes
for k = [-40:10:-10 -5 -1]
    isosurface(xi * 1e2, yi * 1e2, zi * 1e2, I, k)
end

% Use a diverging colormap for better contrast
colormap([1 - hot; hot])

% Show colour scale
colorbar

% Improve figure appearance
box on
grid on

% Label depth axis
zlabel('[cm]')

% Title describing the PSF location and units
title('PSF at (0,0,3) cm [dB]')
