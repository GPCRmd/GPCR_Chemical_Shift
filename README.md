# Integrative Biophysical Illumination of the 3D GPCRome Dynamics

**Authors:**  
Brian Medel-Lacruz, David Aranda-García, Aditya Prasad Patra, Tomasz Maciej Stepniewski, Mariona Torrens-Fontanals, Franz Hagn, Massimiliano Bonomi, Jiafei Mao, Jana Selent

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

---

**Contact:**  
[Brian Medel](mailto:brianmedelmo@gmail.com)
