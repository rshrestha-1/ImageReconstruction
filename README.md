# Image Reconstruction Project

This ImageReconstruction repository contains the Matlab/Python code developed for the data processing and image reconstruction pipelines in 
decoding the integration of the acoustic lens with the L12-3v US linear probe.
The repository is organised into separate packages for various stages and aims of the project.

---

## Repository Structure

```text
project-root/
│── data/
│── preprocessing/
│── reconstruction/
│── models/
│── evaluation/
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

#### `/data`
Description of what is stored here.

Example:
Contains datasets used for training, testing, or reconstruction. This may include raw images, processed images, and sample inputs. Large datasets may be excluded from GitHub and accessed separately.

#### `/preprocessing`
Description.

Example:
Contains scripts for image cleaning, normalisation, resizing, filtering, augmentation, or preparing input data before reconstruction.

#### `/reconstruction`
Description.

Example:
Core reconstruction algorithms and pipelines. Includes methods used to reconstruct images from degraded, incomplete, or transformed inputs.

#### `/models`
Description.

Example:
Stores trained models, model architectures, or saved checkpoints used during reconstruction.

#### `/evaluation`
Description.

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
