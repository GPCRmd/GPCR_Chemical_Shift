# Illuminating underexplored GPCR dynamics via GPCRMD chemical shift streaming  - Neurotensin Receptor

This repository contains the data, scripts, and analysis for the section: **"Illuminating underexplored GPCR dynamics via GPCRMD chemical shift streaming - Neurotensin Receptor."** 

---

## Repository Structure

The repository is organized into the following folders:

### 1. **Computational_Data**
Contains the predicted CS values for the neurotensin receptor complexed with the peptide ligand. The folder includes:

- **Structural files (PDB)**: Structures used for MD simulations.
  - `protein.pdb`

- **Predicted CS values**: CSV files containing the predicted chemical shifts. Download them from [here](https://www.gpcrmd.org/view/2099/) and paste them into this folder to run the notebook (rename them as follows):
  - `cs_rep1.csv`
  - `cs_rep2.csv`
  - `cs_rep3.csv`

---

### 2. **Plots**
Output plots generated with the notebook. They include specific plots for each methionine and the histograms comparing them.

---

### 3. **Notebook.ipynb**
A Jupyter notebook with the workflow used for this section.

---

