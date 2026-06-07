
% Example:KWAVE
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
param.fc = 7.54e6;

% Fractional bandwidth of the transducer in percent
param.bandwidth = 93;

% Width of each transducer element [m]
param.width = 1.7000e-04;

% Height of each transducer element [m]
param.height = 0.0050;

% Sampling frequency of the transducer (30 MHz)
param.fs = 4*param.fc;

% -------------------------------------------------------------------------
% MATRIX ARRAY GEOMETRY
% -------------------------------------------------------------------------

% Pitch (centre-to-centre spacing) between elements [m]
param.pitch = 2.0000e-04;

% Set the number of elements in the transducer array
param.Nelements = 192;

% The grid is centred around (0,0)

% x positions (centred linear array)
xe = ((1:param.Nelements) - (param.Nelements+1)/2) * param.pitch;

% y positions (all zero for linear probe)
ye = zeros(1, param.Nelements);

% Store element positions as a 2 x Nelements matrix:
%   row 1 -> x positions
%   row 2 -> y positions
param.elements = [xe; ye];

% -------------------------------------------------------------------------
% POINT SCATTERER (PSF) LOCATION
% -------------------------------------------------------------------------

% Position of the point scatterer [m]
x0 = 0;        % lateral x
y0 = 0;        % lateral y
z0 = 4e-2;     % depth (4 cm)

% -------------------------------------------------------------------------
% TRANSMIT SETTINGS
% -------------------------------------------------------------------------

% Transmit time delays for plane-wave transmission
% (all zeros = plane wave)
dels = zeros(1, param.Nelements);   % N elements

% -------------------------------------------------------------------------
% RF DATA SIMULATION
% -------------------------------------------------------------------------

RF = zeros(305, param.Nelements);  % base of empty RF data

% Replace elements of choice into the RF
RF(:,87) = sensor_data_filtered(:,1);
RF(:,97) = sensor_data_filtered(:,2);
RF(:,106) = sensor_data_filtered(:,3);
%% 

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
xi = voxel_map_x;
yi = voxel_map_y;
zi = voxel_map_z;

% -------------------------------------------------------------------------
% 3D DELAY-AND-SUM BEAMFORMING
% -------------------------------------------------------------------------

% Perform 3D DAS beamforming:
%   - applies geometric delays
%   - sums contributions from all elements
% Result is a complex-valued 3D image
IQb = das(IQ,xi,zi,dels,param);

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
