# Integrative Biophysical Illumination of the 3D GPCRome Dynamics

**Authors:**  
Brian Medel-Lacruz, Adrián García-Recio, David Aranda-García, Aditya Prasad Patra, Tomasz Maciej Stepniewski, Mariona Torrens-Fontanals, Franz Hagn, Massimiliano Bonomi, Jiafei Mao, Jana Selent

## Repository Structure

In the following folders, you will find the necessary scripts and data to reproduce the results and figures described in the publication.

Due to size constraints, **molecular dynamics (MD) trajectories** and some **chemical shift (CS) prediction files** are not included in this repository. However, you can download these files directly from [GPCRmd](https://gpcrmd.org/).

Each folder contains:
- A `README.md` with detailed explanations of how to download the data
- A `Jupyter Notebook` to guide you step-by-step through the workflows

## Environment Setup

To reproduce the analyses in this repository, please use the provided Conda environment file:

```bash
conda env create -f CS_environment.yml
conda activate CS_env
```

*Note: Make sure to have Conda installed before running these commands.*

# Chemical Shift prediction (CS_computation)
*Note: the chemical shift prediction scripts were built for our system machine that contains specific rutes and it is necessary to modify them to include input and output paths.

The chemical shift prediction is performed using two ways:
- Mdtraj: SPARTA+ & Shiftyx2
- UCBshift: CSpred  

## Environment Setup
To run both ways is necessary to install the next software: 
- DSSP - v.2.04
- BLAST+ - v.2.9.0
- mtm-align - v.20180725
- reduce - v.3.23
- CSpred: becomes from UCBShift2.0 program: https://github.com/THGLab/CSpred 
- Miniconda3

To import the environment a .yml file is generated 'CS_computation/cscomp_environment.yml'

```bash
conda env create -f cscomp_environment.yml
conda activate cspred
```

## Run program

On both ways is necessary to configure the path of input and output files, and the software program paths. 
### SPARTA+ & SHIFTYX2

```bash
nohup python get_chemshift.py --dynid 2375 --cores 3 > cs_sparta_shifty.log 2>&1 & 
```

### UCBSHIFT

```bash
nohup python get_cspred.py --dynid 2375 --cores 3 > cs_ucbshift.log 2>&1 & 
```

## Example input
To run the scripts, we used the following files stored in GPCRmd: https://www.gpcrmd.org/dynadb/dynamics/id/2375/ and the compl_info.json file.

*Note: The trajectory .dcd file have a big size. For this reason, it is necessary to use the link to download it. 

## Example output

The directory example contains output files from these scripts used: 

- SPARTA+: cs_sparta_* 
- Shiftyx2: cs_shifty_*
- UCBshift: cs_ucbshift_*

---

**Contact:**  
[Brian Medel](mailto:brianmedelmo@gmail.com)
