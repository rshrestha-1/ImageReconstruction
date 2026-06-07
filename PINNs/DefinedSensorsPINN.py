## Imports
# Deep learning framework
import torch
import torch.nn as nn

# Numerical handling
import numpy as np

# MATLAB file handling
from scipy.io import loadmat, savemat

import random

# Set seeds for reproducibility
seed = 42
torch.manual_seed(seed)            # PyTorch CPU
torch.cuda.manual_seed(seed)       # PyTorch GPU (if using)
torch.cuda.manual_seed_all(seed)   # All GPUs
np.random.seed(seed)               # NumPy
random.seed(seed)                  # Python built-in RNG

## Load
# Load your MATLAB file
data = loadmat('/Users/rebeccashrestha/Library/CloudStorage/OneDrive-ImperialCollegeLondon/DAPP/Y3/Matlab/Sensor data/H_exp_new.mat')

# Original H shape: (475136, 216)
H = data['H_exp_new']

M, N = H.shape

# Known dataset structure
num_sensors = 58
num_time = 8192
Nx = Ny = Nz = 6

print("Loaded H:", H.shape)

## Reshape H to tensor
# Convert from (475136, 216) → (58, 8192, 216)
H_reshaped = H.reshape(num_sensors, num_time, N)

print("Reshaped H:", H_reshaped.shape)

## FFT over time dimension

dt = 1/30e6  # Sampling interval (s)(30 MHz sampling rate)

# Perform FFT along time axis
H_fft = np.fft.fft(H_reshaped, axis=1)

# Frequency values
freqs = np.fft.fftfreq(num_time, d=dt)

# Keep only positive frequencies
mask = freqs > 0
H_fft = H_fft[:, mask, :]
freqs = freqs[mask]

print("FFT done. Shape:", H_fft.shape)

## Select frequencies of interest (based on energy or known frequencies)
# Compute energy per frequency
energy = np.sum(np.abs(H_fft)**2, axis=(0,2))

# Select top K frequencies
K = 20
top_idx = np.argsort(energy)[-K:]

H_fft = H_fft[:, top_idx, :]
freqs = freqs[top_idx]

print("Selected frequencies:", freqs)

## Build voxel coordinates for input space
coords = []

for idx in range(N):
    z = idx // (Nx * Ny)
    rem = idx % (Nx * Ny)
    y = rem // Nx
    x = rem % Nx
    coords.append([x, y, z])

coords = np.array(coords, dtype=np.float32)

# Normalise to [-1,1]
coords = 2 * (coords - coords.min(0)) / (coords.max(0) - coords.min(0)) - 1

coords = torch.tensor(coords)

## Sensor and frequency normalisation for input space
# Sensor indices
sensor_vals_raw = np.concatenate([
    np.arange(1, 21),     # 1 to 20
    np.arange(46, 61),    # 46 to 60
    np.arange(85, 108)    # 85 to 107
], axis=0)

assert len(sensor_vals_raw) == 58, "Expected 58 sensor indices, received {}".format(len(sensor_vals_raw))

# Normalise to [-1, 1]
sensor_vals = 2 * (sensor_vals_raw - sensor_vals_raw.min()) / (sensor_vals_raw.max() - sensor_vals_raw.min()) - 1

# Convert to torch tensor
sensor_vals = torch.tensor(sensor_vals, dtype=torch.float32)

# Normalise frequencies
f_norm = 2 * (freqs - freqs.min()) / (freqs.max() - freqs.min()) - 1
f_norm = torch.tensor(f_norm, dtype=torch.float32)

## Wave number k(f) calculation
# Define wave speed
c = 1500  # (m/s)

k_vals = 2 * np.pi * freqs / c
k_vals = torch.tensor(k_vals, dtype=torch.float32)

## Define Physics-Informed Neural Network (PINN)
class PINN(nn.Module):
    def __init__(self):
        super().__init__()

        n_fourier = 128  # number of Fourier features

        # Random projection matrix (fixed, not trained)
        self.register_buffer('B', torch.randn(3, n_fourier) * 5.0)
        # scaling (5.0) controls frequency range — can tune later

        input_dim = 2 * n_fourier + 2  # sin + cos + s + f

        # Hidden layers
        self.net = nn.Sequential(
            nn.Linear(input_dim, 128),
            nn.Tanh(),

            nn.Linear(128, 128),
            nn.Tanh(),

            nn.Linear(128, 128),
            nn.Tanh(),

            nn.Linear(128, 2)  # real + imag
        )

    def forward(self, x, y, z, s, f):

        # Stack spatial coords
        xyz = torch.stack([x, y, z], dim=1)  # (N, 3)

        # Fourier feature projection
        proj = xyz @ self.B  # (N, n_fourier)

        # Apply sin and cos
        ffe = torch.cat([torch.sin(proj), torch.cos(proj)], dim=1)

        # Concatenate with sensor and frequency (NOT encoded)
        inputs = torch.cat([ffe, s.unsqueeze(1), f.unsqueeze(1)], dim=1)

        # Feed through network
        out = self.net(inputs)

        return out[:,0], out[:,1]
    
model = PINN()

## Derivative function
def grad(u, x):
    return torch.autograd.grad(
        u, x,
        grad_outputs=torch.ones_like(u),
        create_graph=True
    )[0]

## Physics loss function
def physics_loss(model, x, y, z, s, f, k):

    x.requires_grad_(True)
    y.requires_grad_(True)
    z.requires_grad_(True)

    u_r, u_i = model(x, y, z, s, f)

    # First derivatives
    u_r_x = grad(u_r, x)
    u_r_y = grad(u_r, y)
    u_r_z = grad(u_r, z)

    u_i_x = grad(u_i, x)
    u_i_y = grad(u_i, y)
    u_i_z = grad(u_i, z)

    # Second derivatives
    u_r_xx = grad(u_r_x, x)
    u_r_yy = grad(u_r_y, y)
    u_r_zz = grad(u_r_z, z)

    u_i_xx = grad(u_i_x, x)
    u_i_yy = grad(u_i_y, y)
    u_i_zz = grad(u_i_z, z)

    # Helmholtz residual
    res_r = u_r_xx + u_r_yy + u_r_zz + (k**2) * u_r
    res_i = u_i_xx + u_i_yy + u_i_zz + (k**2) * u_i

    return torch.mean(res_r**2 + res_i**2)

## Training loop
optimizer = torch.optim.Adam(model.parameters(), lr=1e-3)

batch_size = 8000
loss_total_hist = []
loss_data_hist = []
loss_phys_hist = []

for epoch in range(3000):

    optimizer.zero_grad()

    # Random sampling
    s_idx = torch.randint(0, len(sensor_vals), (batch_size,))
    f_idx = torch.randint(0, K, (batch_size,))
    v_idx = torch.randint(0, N, (batch_size,))

    # Inputs
    xyz = coords[v_idx]
    x, y, z = xyz[:,0], xyz[:,1], xyz[:,2]

    s = sensor_vals[s_idx]
    f = f_norm[f_idx]
    k = k_vals[f_idx]

    # Ground truth
    H_sample = H_fft[s_idx, f_idx, v_idx]
    H_r = torch.tensor(np.real(H_sample), dtype=torch.float32)
    H_i = torch.tensor(np.imag(H_sample), dtype=torch.float32)

    # Data loss
    u_r, u_i = model(x, y, z, s, f)
    loss_d = torch.mean((u_r - H_r)**2 + (u_i - H_i)**2)

    # Physics loss
    loss_p = physics_loss(model, x, y, z, s, f, k)

    loss = loss_d + 1e-3 * loss_p

    loss_total_hist.append(loss.item())
    loss_data_hist.append(loss_d.item())
    loss_phys_hist.append(loss_p.item())

    loss.backward()
    optimizer.step()

    if epoch % 200 == 0:
        print(f"Epoch {epoch}, Loss: {loss.item():.6f}")

## Generate new voxel grid and predict with trained PINN
# Original physical size
phys_size_orig = 6.0  # mm
Nx_orig = Ny_orig = Nz_orig = 6

# New physical size (extrapolation region)
phys_size_new = 10.0  # mm
Nx_new = Ny_new = Nz_new = 20  # Increase voxel resolution for smoother prediction

# Generate new voxel coordinates in physical units
coords_new = []

for z in range(Nz_new):
    for y in range(Ny_new):
        for x in range(Nx_new):
            coords_new.append([x, y, z])

coords_new = np.array(coords_new, dtype=np.float32)

# Normalise coordinates to [-1,1] based on the new physical size
coords_new[:,0] = 2 * (coords_new[:,0] / (Nx_new - 1) * phys_size_new / phys_size_orig) - 1
coords_new[:,1] = 2 * (coords_new[:,1] / (Ny_new - 1) * phys_size_new / phys_size_orig) - 1
coords_new[:,2] = 2 * (coords_new[:,2] / (Nz_new - 1) * phys_size_new / phys_size_orig) - 1

coords_new = torch.tensor(coords_new, dtype=torch.float32)

## Build new H
H_new = np.zeros((num_sensors, K, len(coords_new)), dtype=np.complex64)

for s_idx in range(num_sensors):

    s_val = sensor_vals[s_idx]

    for f_idx in range(K):

        f_val = f_norm[f_idx]

        s_vec = torch.ones(len(coords_new)) * s_val
        f_vec = torch.ones(len(coords_new)) * f_val

        u_r, u_i = model(coords_new[:,0], coords_new[:,1], coords_new[:,2], s_vec, f_vec)

        H_new[s_idx, f_idx, :] = u_r.detach().numpy() + 1j*u_i.detach().numpy()

## Rebuild frequency domain data on new grid

# Total time samples
Nt = 8192

# Create empty full spectrum
# Shape: (sensors, full_freqs, voxels)
H_full_fft = np.zeros((num_sensors, Nt, len(coords_new)), dtype=np.complex64)

# Get original FFT frequency indices
full_freqs = np.fft.fftfreq(Nt, d=dt)

# Find where our selected freqs sit in full spectrum
for i, f in enumerate(freqs):

    # Find closest index in full frequency array
    idx = np.argmin(np.abs(full_freqs - f))

    # Insert learned values
    H_full_fft[:, idx, :] = H_new[:, i, :]

    # Also insert complex conjugate for negative frequencies
    neg_idx = np.argmin(np.abs(full_freqs + f))
    H_full_fft[:, neg_idx, :] = np.conj(H_new[:, i, :])

## Inverse FFT

# Apply inverse FFT along time axis
H_time = np.fft.ifft(H_full_fft, axis=1)

print("Reconstructed time-domain shape:", H_time.shape)

## Reshape to original matrix form

# Currently: (58, 8192, new_voxels)
# Convert back to (475136, new_voxels)
H_time_matrix = H_time.reshape(num_sensors * Nt, len(coords_new))

## Save final time-domain extrapolated H

savemat('H_extrapolated_time.mat', {
    'H_time': H_time_matrix,
    'Nx': Nx_new,
    'Ny': Ny_new,
    'Nz': Nz_new
})

print("Saved final time-domain extrapolated H")