# Integrative enhancement of spectroscopic pipelines via GPCRMD-CS - rhodopsin receptor

This repository contains the data, scripts, and analysis for the section: **"Integrative enhancement of spectroscopic pipelines via GPCRMD-CS - rhodopsin receptor."** 

---

## Repository Structure

The repository is organized into the following folders:

### 1. **Computational_Data**
Contains the predicted CS values for the rhodopsin receptor. The folder includes:

- **Structural files (PDB)**: Structures used for MD simulations.
  - `17954_dyn_872.pdb`

- **Trajectory files (XTC)**: XTC file for the system. Download it from [here](https://www.gpcrmd.org/view/872/) and paste them into this folder to run the notebook.
  - `16264_trj_872.xtc`

- **Predicted CS values**: CSV file containing the predicted chemical shifts. Download it from [here](https://www.gpcrmd.org/view/872/) and paste them into this folder to run the notebook.
  - `cs_dyn872_16264.csv`

---

### 2. **Notebook.ipynb**
A Jupyter notebook with the workflow used for this section.

---