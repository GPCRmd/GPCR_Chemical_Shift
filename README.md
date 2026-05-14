# Integrative Biophysical Illumination of the 3D GPCRome Dynamics

**Authors:**  
Brian Medel-Lacruz, Adrián García-Recio, David Aranda-García, Aditya Prasad Patra, Tomasz Maciej Stepniewski, Mariona Torrens-Fontanals, Franz Hagn, Massimiliano Bonomi, Jiafei Mao, Jana Selent

## Study Overview & Data Statistics:

The computational analysis presented in this study represents a large-scale characterization of the GPCRome. The total data processed includes:

* **Unique Dynamics Analyzed**: 1,028
* **Total Trajectories**: 3,034
* **Total Frames Analyzed**: 4,876,719
* **Cumulative Data Size**: 255 GB

## Repository Structure & Data Access:

**For General Analysis:**
In the following folders, you will find the necessary scripts and data to reproduce the results and figures described in the publication.

Each folder contains:
- A `README.md` with detailed explanations of how to download the data.
- A `Jupyter Notebook` to guide you step-by-step through the workflows.

**For Chemical Shift Computation:**

Due to the large-scale nature of the GPCRome dataset, **molecular dynamics (MD) trajectories** and **full chemical shift (CS) prediction files** are hosted externally on the [GPCRmd platform](https://gpcrmd.org/).

To execute the provided scripts and reproduce the analysis, please follow these guidelines:

* **Required Data**: The analysis requires the compl_info.json file and dynamics data. As a reference example, we utilize Dynamics ID 2375, which can be accessed here: [GPCRmd Dynamics 2375](https://www.gpcrmd.org/dynadb/dynamics/id/2375/).
* **Trajectory Handling**: Note that .dcd trajectory files possess significant file sizes. We recommend downloading these files directly via the GPCRmd links before running the local prediction workflows.

---

## Software & Environment Requirements:

The analysis was performed using a suite of chemical shift prediction and trajectory processing tools integrated within a Miniconda3 environment.
Tested Software Versions

Chemical Shift Predictors:
* **UCBShift v2.0 (CSpred)**: [Source](https://github.com/THGLab/CSpred)
* **SPARTA+**: [Integrated via MDTraj 1.9.4](https://mdtraj.org/1.9.4/api/generated/mdtraj.chemical_shifts_spartaplus.html#mdtraj.chemical_shifts_spartaplus)
* **ShiftyX2**: [Integrated via MDTraj 1.9.4](https://mdtraj.org/1.9.4/api/generated/mdtraj.chemical_shifts_shiftx2.html#mdtraj.chemical_shifts_shiftx2)

Core Dependencies:
* **MDTraj 1.9.4**: For trajectory handling and structural complements.
* **DSSP v2.04**: Secondary structure assignment.
* **BLAST+ v2.9.0**: Sequence alignment.
* **mtm-align v20180725**: Structural alignment.
* **Reduce v3.23**: Hydrogen addition and optimization.

---

## Installation & Benchmarking:

To ensure reproducibility, we provide the following estimates for a "normal" desktop computer (standard modern CPU/RAM configuration).

**1. Typical Installation Time:**

* **Estimated Time**: 20–45 minutes.
* **Details**: This includes setting up the Miniconda environment, downloading pre-compiled binaries via Conda, compiling source-dependent tools (e.g., mtm-align), and configuring BLAST databases.

**2. Expected Run Time (Demo Dataset):**

Execution times were measured using the Python `time` module. Benchmarks are provided for different computational loads:

| Analysis Type | Dataset Size | Total Seconds | Formatted Time |
| :--- | :--- | :--- | :--- |
| **SPARTA+** | Standard Demo | 18,207.6 s | ~5h 03m |
| **UCBShift (Short)** | Minimum Frame Set | 121,866.5 s | ~33.8h |
| **UCBShift (Median)** | Average Frame Set | 173,795.8 s | ~48.3h |
| **UCBShift (Long)** | Maximum Frame Set | 195,773.8 s | ~54.4h |

---

## Environment Setup:

To reproduce the analyses in this repository, please use the provided Conda environment file:

**For General Analysis:**

```bash
conda env create -f CS_environment.yml
conda activate CS_env
```

**For Chemical Shift Computation (CS_computation):**

```bash
conda env create -f cscomp_environment.yml
conda activate cspred
```

***Note**: the chemical shift prediction scripts were built for our system machine that contains specific rutes and it is necessary to modify them to include input and output paths.*

---

## Usage:

***Note**: Prediction scripts include specific system paths. Users must modify the input and output paths in the scripts to match their local environment.*

**Running Predictors**
### SPARTA+ & SHIFTYX2:

```bash
nohup python get_chemshift.py --dynid 2375 --cores 3 > cs_sparta_shifty.log 2>&1 & 
```

### UCBSHIFT(CSpred):

```bash
nohup python get_cspred.py --dynid 2375 --cores 3 > cs_ucbshift.log 2>&1 & 
```

---

**Contact:**  
[Brian Medel](mailto:brianmedelmo@gmail.com)
