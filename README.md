# Image Reconstruction Project

This ImageReconstruction repository contains the Matlab/Python code developed for the data processing and image reconstruction pipelines in decoding the integration of the acoustic lens with the Ultrasound L12-3v linear probe.
The repository is organised into separate packages for various stages and aims of the project.

---

## Repository Structure

```text
project-root/
│── 2d_das/
│── 3D_trials/
│── Algorithms/
│── Execution/
│── Experimental/
│── February/
│── MUST/
│── PINNs/
│── ezdaSS/
│── figures/
│── README.md
│── requirements.txt
```

### Folder Overview

#### `/2d_das`
Contains testing from the 2D DAS reconstruction period early on within autumm, requires installation of MUST modules.

#### `/3D_trials`
Contains testing from 3D DAS and moving into 3D reconstruction territory also early on, with initial trial u, H, v data.

#### `/Algorithms`
Core reconstruction algorithms and descriptions for function use (inputs/outouts). Includes main methods used to reconstruct images.

#### `/Execution`
Scripts used to execute the reconstruction algorithms. MATLAB sections numbered and referenced both between and within scripts. Contains scripts for normalisation, resizing, filtering, autoconvolution, and preparing input data for reconstruction. Contains various unit testing of methods and sanity tests.

#### `/Experimental`
Contains code for processing the limited experimental data. Includes ASA.

Example:
Scripts for assessing reconstruction quality using metrics such as PSNR, SSIM, MSE, or qualitative comparisons.

#### `/February`
Contains scripts to debug and run the original code written up prior to receiving data, once the first sets of simulation data was being provided.

Example:
Contains plotting scripts and tools for displaying reconstructed images, comparisons, graphs, and performance metrics.

#### `/MUST`
Helper functions from MUST shared across the repository.

#### `/PINNs`
Contains the Python scripts to run PINNs to extrapolate/interpolate the limited experimental data. Includes versions with/without ffe; extraction with .mat/HDF5; reduced hidden layers/adaptive physics loss weight etc. to reduce instability burdens. Stores model architectures.

#### `/ezdaSS`
Contains further MUST modules.

#### `/figures`
Contains figure visualisation code.

---

## Dependencies

- Python `x.x`
- NumPy
- SciPy
- PyTorch

Install all dependencies using:

```bash
pip install -r requirements.txt
```
---

## License

This project is intended for academic use.
