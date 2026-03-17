import os
import pandas as pd
import numpy as np
import statistics
import matplotlib.pyplot as plt
from scipy.stats import pearsonr

# === 1. Create output directory ===
output_dir = "results_alltypes"
os.makedirs(output_dir, exist_ok=True)
print(f"All results will be saved in: {output_dir}/")

# === Functions ===
def read_experimental_cs(file_path, shift=0, resid_column=19, cs_column=10, atom_column=7, resname_column=6):
    data_dict = {}
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split()
            resname = columns[resname_column]
            resid = int(columns[resid_column]) + shift
            atomname = columns[atom_column]
            chemshift = float(columns[cs_column])
            atomname_base = ''.join(filter(str.isalpha, atomname))
            data_dict[f"{resid}_{resname}_{atomname_base}"] = chemshift
    return data_dict

def extract_residue_ids(exp_dict):
    return {int(k.split('_')[0]) for k in exp_dict.keys()}

def read_computational_cs(file_paths, residue_ids):
    dict_comp_values = {}
    for file_path in file_paths:
        with open(file_path, 'r') as file:
            lines = file.readlines()[1:]
            for line in lines:
                data = line.strip().split(";")
                if int(data[0]) not in residue_ids:
                    continue
                float_list = [float(element) for element in data[4:-1]]
                atomname = data[1]
                atomname_base = ''.join(filter(str.isalpha, atomname))
                id_key = f"{int(data[0])}_{data[2]}_{atomname_base}"
                dict_comp_values[id_key] = round(statistics.mean(float_list), 5)
    return dict_comp_values

def combine_replicates_comp(list_of_dicts):
    combined_dict = {}
    for d in list_of_dicts:
        for k, v in d.items():
            combined_dict.setdefault(k, []).append(v)
    return {k: round(statistics.mean(v), 5) for k, v in combined_dict.items()}

def read_81_proves_experimental(file_path):
    data_dict = {}
    with open(file_path, 'r') as file:
        for line in file:
            columns = line.strip().split()
            if len(columns) < 2:
                continue
            parts = columns[0].split('.')
            if len(parts) < 2:
                continue
            resid = parts[0]
            atomname = parts[-1]
            chemshift = float(columns[1])
            atomname_base = ''.join(filter(str.isalpha, atomname))
            if atomname_base in ['H', 'N']:
                data_dict[f"{resid}_{atomname_base}"] = chemshift
    return data_dict

def convert_index_comp_dict_81(comp_dict):
    converted = {}
    for key, val in comp_dict.items():
        parts = key.split('_')
        if len(parts) == 3:
            resid, _, atomname = parts
            atom_base = ''.join(filter(str.isalpha, atomname))
            if atom_base in ['H', 'N']:
                converted[f"{resid}_{atom_base}"] = val
    return converted

# === Plot functions ===
def plot_by_atomtype(list_exp_dict, list_comp_dict, title, output_path, atom_type, exact_match=False):
    exp_shifts, comp_shifts = [], []
    for exp_dict, comp_dict in zip(list_exp_dict, list_comp_dict):
        for key in exp_dict.keys():
            if key in comp_dict:
                atom_name = key.split('_')[-1]
                if (exact_match and atom_name == atom_type) or (not exact_match and atom_name.startswith(atom_type)):
                    exp_shifts.append(exp_dict[key])
                    comp_shifts.append(comp_dict[key])

    if len(exp_shifts) < 2:
        return None

    exp_arr, comp_arr = np.array(exp_shifts), np.array(comp_shifts)
    # Compute Pearson correlation properly
    r, _ = pearsonr(comp_arr, exp_arr)
    rmse = np.sqrt(np.mean((exp_arr - comp_arr) ** 2))

    plt.figure(figsize=(4, 4))
    plt.scatter(comp_arr, exp_arr, alpha=0.6, color='black')
    plt.plot([min(comp_arr), max(comp_arr)], [min(comp_arr), max(comp_arr)], 'r--')
    plt.xlabel(f'Computational CS (ppm)', size=14)
    plt.ylabel(f'Experimental CS (ppm)', size=14)
    plt.title(title)
    plt.text(0.05, 0.95, f"r = {r:.3f}\nRMSE = {rmse:.3f}",
             transform=plt.gca().transAxes, fontsize=11, va='top',
             bbox=dict(facecolor='white', alpha=0.8))
    plt.tight_layout()
    plt.grid(True)
    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)
    plt.savefig(output_path)
    plt.close()

    return {"atom_type": atom_type, "pearson_r": r, "rmse": rmse, "n_points": len(exp_shifts)}

def plot_by_atomtype_color_each_experiment(list_exp_dict, list_comp_dict, title, output_path, atom_type, exact_match=False):
    colors = ['red', 'blue', 'green', 'orange', 'purple', 'brown', 'pink', 'gray', 'olive', 'cyan']
    plt.figure(figsize=(4, 4))
    all_exp, all_comp = [], []

    for i, (exp_dict, comp_dict) in enumerate(zip(list_exp_dict, list_comp_dict)):
        exp_shifts, comp_shifts = [], []
        for key in exp_dict.keys():
            if key in comp_dict:
                atom_name = key.split('_')[-1]
                if (exact_match and atom_name == atom_type) or (not exact_match and atom_name.startswith(atom_type)):
                    exp_shifts.append(exp_dict[key])
                    comp_shifts.append(comp_dict[key])
        if exp_shifts:
            plt.scatter(comp_shifts, exp_shifts, alpha=0.6, color=colors[i % len(colors)], label=f'Exp {i+1}')
            all_exp.extend(exp_shifts)
            all_comp.extend(comp_shifts)

    if len(all_exp) > 1:
        exp_arr, comp_arr = np.array(all_exp), np.array(all_comp)
        r, _ = pearsonr(comp_arr, exp_arr)
        rmse = np.sqrt(np.mean((exp_arr - comp_arr) ** 2))
        plt.plot([min(comp_arr), max(comp_arr)], [min(comp_arr), max(comp_arr)], 'r--', label='y=x')
        plt.text(0.05, 0.95, f"r = {r:.3f}\nRMSE = {rmse:.3f}",
                 transform=plt.gca().transAxes, fontsize=11, va='top',
                 bbox=dict(facecolor='white', alpha=0.8))

    plt.xlabel(f'Computational CS (ppm)', size=14)
    plt.ylabel(f'Experimental CS (ppm)', size=14)
    plt.title(title)
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.xticks(fontsize=11)
    plt.yticks(fontsize=11)
    plt.savefig(output_path)
    plt.close()

# === Load experimental data ===
data_dict_exp1 = read_experimental_cs("systems/exp1/cs_data.txt")
sel_resid1 = extract_residue_ids(data_dict_exp1)
data_dict_comp1 = read_computational_cs(["systems/exp1/cs_dyn143_11356.csv"], sel_resid1)

data_dict_exp2 = read_81_proves_experimental("systems/exp2/shifts_81pro.txt")
sel_resid2 = extract_residue_ids(data_dict_exp2)
data_dict_comp2 = combine_replicates_comp([
    read_computational_cs(["systems/exp2/cs_dyn195_11792.csv"], sel_resid2),
    read_computational_cs(["systems/exp2/cs_dyn195_11793.csv"], sel_resid2),
    read_computational_cs(["systems/exp2/cs_dyn195_11846.csv"], sel_resid2)
])
data_dict_comp2_converted = convert_index_comp_dict_81(data_dict_comp2)

data_dict_exp3 = read_experimental_cs("systems/exp3/cs_data.txt", -9)
sel_resid3 = extract_residue_ids(data_dict_exp3)
data_dict_comp3 = combine_replicates_comp([
    read_computational_cs(["systems/exp3/cs_dyn44_10462.csv"], sel_resid3),
    read_computational_cs(["systems/exp3/cs_dyn44_10463.csv"], sel_resid3),
    read_computational_cs(["systems/exp3/cs_dyn44_10464.csv"], sel_resid3)
])

data_dict_exp5 = read_experimental_cs("systems/exp5/cs_data.txt", 0, 20, 11, 8, 7)
sel_resid5 = extract_residue_ids(data_dict_exp5)
data_dict_comp5 = combine_replicates_comp([
    read_computational_cs(["systems/exp5/cs_dyn2359_25636.csv"], sel_resid5),
    read_computational_cs(["systems/exp5/cs_dyn2359_25637.csv"], sel_resid5),
    read_computational_cs(["systems/exp5/cs_dyn2359_25638.csv"], sel_resid5)
])

data_dict_exp6 = read_experimental_cs("systems/exp6/cs_data.txt", 0, 20, 11, 8, 7)
sel_resid6 = extract_residue_ids(data_dict_exp6)
rep6_files = [f"systems/exp6/cs_dyn1753_223{i}.csv" for i in range(37, 47)]
data_dict_comp6 = combine_replicates_comp([read_computational_cs([f], sel_resid6) for f in rep6_files])

# === Combine all experiments ===
list_exp_dict = [data_dict_exp1, data_dict_exp2, data_dict_exp3, data_dict_exp5, data_dict_exp6]
list_comp_dict = [data_dict_comp1, data_dict_comp2_converted, data_dict_comp3, data_dict_comp5, data_dict_comp6]

# === Compute & plot correlations by atom type ===
atom_types_exact = ["C", "CA", "CB", "N", "H"]
results = []

for atom in atom_types_exact:
    out_path = os.path.join(output_dir, f"all_experiments_{atom}.png")
    res = plot_by_atomtype(list_exp_dict, list_comp_dict, f"All Experiments: {atom} atoms",
                           out_path, atom, exact_match=True)
    if res:
        results.append(res)

    out_path_color = os.path.join(output_dir, f"all_experiments_{atom}_color.png")
    plot_by_atomtype_color_each_experiment(list_exp_dict, list_comp_dict,
                                           f"All Experiments (Colored): {atom} atoms",
                                           out_path_color, atom, exact_match=True)

# === Compute and plot per-element correlations ===
element_groups = {
    "C": ["C", "CA", "CB", "CE", "CD"],
    "N": ["N"],
    "H": ["H"]
}

element_results = []
for elem, subtypes in element_groups.items():
    all_exp, all_comp = [], []
    color_exp_data = []

    for i, (exp_dict, comp_dict) in enumerate(zip(list_exp_dict, list_comp_dict)):
        exp_e, comp_e = [], []
        for key in exp_dict.keys():
            if key in comp_dict:
                atomtype = key.split('_')[-1]
                if atomtype in subtypes:
                    exp_e.append(exp_dict[key])
                    comp_e.append(comp_dict[key])
        if exp_e:
            color_exp_data.append((i, exp_e, comp_e))
            all_exp.extend(exp_e)
            all_comp.extend(comp_e)

    if len(all_exp) > 1:
        exp_arr, comp_arr = np.array(all_exp), np.array(all_comp)
        r, _ = pearsonr(comp_arr, exp_arr)
        rmse = np.sqrt(np.mean((exp_arr - comp_arr) ** 2))
        element_results.append({"element": elem, "pearson_r": r, "rmse": rmse, "n_points": len(all_exp)})
        
        # if the correlation is rounded to 1.000, set it to 0.999 to avoid plotting issues
        if round(r, 3) == 1.000:
            r = 0.999

        # === Plot (combined) ===
        plt.figure(figsize=(4, 4))
        plt.scatter(comp_arr, exp_arr, alpha=0.6, color='black')
        plt.plot([min(comp_arr), max(comp_arr)], [min(comp_arr), max(comp_arr)], 'r--')
        plt.title(f"All Experiments: {elem} Elements")
        plt.xlabel(f'Computational CS (ppm)', size=14)
        plt.ylabel(f'Experimental CS (ppm)', size=14)
        plt.text(0.05, 0.95, f"r = {r:.3f}\nRMSE = {rmse:.3f}",
                 transform=plt.gca().transAxes, fontsize=11, va='top',
                 bbox=dict(facecolor='white', alpha=0.8))
        plt.tight_layout()
        plt.grid(True)
        plt.xticks(fontsize=11)
        plt.yticks(fontsize=11)
        plt.savefig(os.path.join(output_dir, f"all_experiments_{elem}_element.png"))
        plt.close()

        # === Plot (colored by experiment) ===
        plt.figure(figsize=(4, 4))
        colors = ['red', 'blue', 'green', 'orange', 'purple']
        for i, exp_e, comp_e in color_exp_data:
            plt.scatter(comp_e, exp_e, alpha=0.6, color=colors[i % len(colors)], label=f'Exp {i+1}')
        plt.plot([min(comp_arr), max(comp_arr)], [min(comp_arr), max(comp_arr)], 'r--', label='y=x')
        plt.title(f"All Experiments (Colored): {elem} Elements")
        plt.xlabel(f'Computational CS (ppm)', size=14)
        plt.ylabel(f'Experimental CS (ppm)', size=14)
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        plt.xticks(fontsize=11)
        plt.yticks(fontsize=11)
        plt.savefig(os.path.join(output_dir, f"all_experiments_{elem}_element_color.png"))
        plt.close()

# === Save summary tables ===
df_atoms = pd.DataFrame(results)
df_elements = pd.DataFrame(element_results)
df_atoms.to_csv(os.path.join(output_dir, "correlation_per_atomtype.csv"), index=False)
df_elements.to_csv(os.path.join(output_dir, "correlation_per_element.csv"), index=False)

# === Combined correlation matrix across C, N, H families ===
if not df_elements.empty:
    corr_matrix = df_elements.pivot_table(index="element", values="pearson_r")
    corr_matrix.to_csv(os.path.join(output_dir, "correlation_family_summary.csv"))
    print("\nSaved correlation matrix for families (C, N, H)")

print("\nSaved results:")
print(f"- {output_dir}/correlation_per_atomtype.csv")
print(f"- {output_dir}/correlation_per_element.csv")
print(f"- {output_dir}/correlation_family_summary.csv (C, N, H family summary)")
print("All plots (standard and colored) are stored in the results/ folder.")
