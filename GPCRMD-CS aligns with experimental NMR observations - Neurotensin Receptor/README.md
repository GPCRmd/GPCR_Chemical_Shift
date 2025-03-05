# GPCRMD-CS aligns with experimental NMR observations - Neurotensin Receptor

This repository contains the data, scripts, and analysis for the section: **"GPCRMD-CS aligns with experimental NMR observations - Neurotensin Receptor."** 

---

## Repository Structure

The repository is organized into the following folders:

### 1. **Computational_Data**
Contains the predicted CS values for the neurotensin receptor complexed with the peptide ligand. The folder includes:

- **Structural files (PDB)**: Structures used for MD simulations.
  - `10661_dyn_66.pdb`

- **Predicted CS values**: CSV files containing the predicted chemical shifts. Download them from [here](https://www.gpcrmd.org/view/2097/) and paste them into this folder to run the notebook.
  - `cs_dyn1753_22337.csv`
  - `cs_dyn1753_22338.csv`
  - `cs_dyn1753_22339.csv`
  - `cs_dyn1753_22340.csv`
  - `cs_dyn1753_22341.csv`
  - `cs_dyn1753_22342.csv`
  - `cs_dyn1753_22343.csv`
  - `cs_dyn1753_22344.csv`
  - `cs_dyn1753_22345.csv`
  - `cs_dyn1753_22346.csv`

---

### 2. **Experimental_Data**
Contains the raw experimental data for the validation of predicted CS values:
- `51907_simulated_hsqc_backbone.csv`
- `shift_neurotensin.txt`

---

### 3. **Plots**
Output plots generated with the notebook:
- `Correlation_C.png`
- `Correlation_CA.png`
- `Correlation_CB.png`
- `correlation_coefficients_table.png`
- `Correlation_H.png`
- `Correlation_N.png`

---

### 4. **Notebook.ipynb**
A Jupyter notebook with the workflow used for this section.

---

