# Integrative enhancement of spectroscopic pipelines via GPCRMD-CS - μ-opioid receptor

This repository contains the data, scripts, and analysis for the section: **"Integrative enhancement of spectroscopic pipelines via GPCRMD-CS - μ-opioid receptor."** 

---

## Repository Structure

The repository is organized into the following folders:

### 1. **Computational_Data**
Contains the predicted CS values for the neurotensin receptor complexed with the peptide ligand. The folder includes:

- **Structural files (PDB)**: Structures used for MD simulations.
  - `11880_dyn_199.pdb`

- **Trajectory files (XTC)**: XTC file for the system. Download it from [here](https://www.gpcrmd.org/view/199/) and paste them into this folder to run the notebook.
  - `11888_trj_199.xtc`

- **Predicted CS values**: CSV file containing the predicted chemical shifts. Download it from [here](https://www.gpcrmd.org/view/199/) and paste them into this folder to run the notebook.
  - `cs_dyn199_11888.csv`

---

### 2. **Reference_Systems**
Contains the two reference PDBs (active/inactive).

---

### 3. **Output**
PDB files with the beta values for comparing residues.

---

### 4. **Notebook.ipynb**
A Jupyter notebook with the workflow used for this section.

---

