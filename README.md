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
│── visualisation/
│── utils/
│── notebooks/
│── results/
│── docs/
│── main.py
│── requirements.txt
│── README.md
```

### Folder Overview

#### `/2d_das`
Contains testing from the 2D DAS reconstruction period early on within autumm, requires installation of MUST modules.

#### `/3D_trials`
Contains testing from 3D DAS and moving into 3D reconstruction territory also early on, with initial trial u, H, v data.

Example:
Contains scripts for image cleaning, normalisation, resizing, filtering, augmentation, or preparing input data before reconstruction.

#### `/Algorithms`
Core reconstruction algorithms and descriptions for function use (inputs/outouts). Includes main methods used to reconstruct images.

#### `/Execution`
Scripts used to execute the reconstruction algorithms. MATLAB sections numbered and referenced both between and within scripts. Contains various unit testing of methods and sanity tests also.

Example:
Stores trained models, model architectures, or saved checkpoints used during reconstruction.

#### `/Experimental`
Contains code for processing the limited experimental data. Includes ASA.

Example:
Scripts for assessing reconstruction quality using metrics such as PSNR, SSIM, MSE, or qualitative comparisons.

#### `/visualisation`
Description.

Example:
Contains plotting scripts and tools for displaying reconstructed images, comparisons, graphs, and performance metrics.

#### `/utils`
Description.

Example:
Helper functions and reusable utilities shared across the repository (e.g., file loading, image operations, common processing functions).

#### `/notebooks`
Description.

Example:
Jupyter notebooks used for experimentation, testing ideas, and exploratory analysis.

#### `/results`
Description.

Example:
Stores output images, logs, figures, reconstruction results, and experiment outputs.

#### `/docs`
Description.

Example:
Project documentation, reports, diagrams, or supplementary explanations.

---

## Dependencies

- Python `x.x`
- NumPy
- OpenCV
- Matplotlib
- SciPy
- PyTorch / TensorFlow (if applicable)

Install all dependencies using:

```bash
pip install -r requirements.txt
```
---

## License

This project is intended for academic use.
