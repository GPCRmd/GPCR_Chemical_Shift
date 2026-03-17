# GPCRMD-CS aligns with experimental NMR observations - Bradykinin 1A receptor

This repository contains the data, scripts, and analysis for the section: **"GPCRMD-CS aligns with experimental NMR observations - Bradykinin 1A receptor."** 

---

## Repository Structure

The repository is organized into the following folders:

### 1. **Computational_data**
Contains MD simulation outputs and structural files for the B1R receptor complexed with the peptide ligand. The folder includes:

- **MD trajectories (XTC)**: 3 replicates of 500 ns simulations. Please, download them from gpcrmd (too large for github), links provided below:
  - [output_rep1.xtc](https://www.gpcrmd.org/view/2095/)
  - [output_rep2.xtc](https://www.gpcrmd.org/view/2095/)
  - [output_rep3.xtc](https://www.gpcrmd.org/view/2095/)

- **Structural files (PDB)**: Structures used for MD simulations.
  - `structure_rep1.pdb`
  - `structure_rep2.pdb`
  - `structure_rep3.pdb`

- **Predicted CS values**: CSV files containing the predicted chemical shifts for the peptide ligand.
  - `CS_peptide_7EIB_1.csv`
  - `CS_peptide_7EIB_2.csv`
  - `CS_peptide_7EIB_3.csv`

---

### 2. **Experimental_data**
Contains the raw experimental data for the validation of predicted CS values:
- `BradReceptor_20141126.txt`

---

### 3. **Plots**
Includes visualizations of the DQ-SQ spectra:
- `DQ-SQ_pro8_replicate1.png`
- `DQ-SQ_pro8_replicate2.png`
- `DQ-SQ_pro8_replicate3.png`

---

### 4. **Notebook.ipynb**
A Jupyter notebook with the workflow used to analyze and extract the images for this experiment.

---
