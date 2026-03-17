import mdtraj as md
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import matplotlib as mpl
from scipy.stats import pearsonr

# -----------------------
# Input
# -----------------------
systems = [
    ("MDs/tmp_dyn_0_2716.pdb", "MDs/25636_trj_2359.xtc"),
    ("MDs/tmp_dyn_0_2716.pdb", "MDs/25637_trj_2359.xtc"),
    ("MDs/tmp_dyn_0_2716.pdb", "MDs/25638_trj_2359.xtc"),
]

CS_data = [
    "CS_data/cs_dyn2359_25636.csv",
    "CS_data/cs_dyn2359_25637.csv",
    "CS_data/cs_dyn2359_25638.csv"
]

res_start, res_end = 349, 390
outdir = "output"
os.makedirs(outdir, exist_ok=True)

# -----------------------
# DSSP mapping
# -----------------------
ss_map = {
    'H': 0, 'E': 1, 'G': 2, 'I': 3,
    'B': 4, 'T': 5, 'S': 6, 'C': 7
}

# -----------------------
# DSSP → Coil extraction
# -----------------------
def compute_coil_vectors(pdb, xtc):
    traj = md.load(xtc, top=pdb)
    dssp = md.compute_dssp(traj, simplified=False)

    res_idx = [
        res.index for res in traj.topology.residues
        if res_start <= res.resSeq <= res_end
    ]

    ss_matrix = dssp[:, res_idx]

    ss_numeric = np.array([
        [ss_map.get(code, 7) for code in frame]
        for frame in ss_matrix
    ])

    coil_dict = {}
    for i in range(ss_numeric.shape[1]):
        residue_number = res_start + i
        coil_dict[residue_number] = (ss_numeric[:, i] == ss_map['C']).astype(int)

    return coil_dict


def read_cs_csv(path):
    return pd.read_csv(path, sep=";")


def extract_region_cs(df_cs, atom_types):
    df_region = df_cs[
        (df_cs["resSeq"] >= res_start) &
        (df_cs["resSeq"] <= res_end) &
        (df_cs["name"].isin(atom_types))
    ].sort_values(["resSeq", "name"])

    frame_cols = [
        c for c in df_region.columns
        if c not in ["resSeq", "name", "resname", "resname_s"]
    ]

    return df_region, frame_cols


# -----------------------
# Compute coil data
# -----------------------
dict_coil_per_residue = {}

for i, (pdb, xtc) in enumerate(systems, start=1):
    label = f"rep{i}"
    dict_coil_per_residue[label] = compute_coil_vectors(pdb, xtc)

# -----------------------
# Concatenate coil
# -----------------------
concatenated_coil = {}

for res in range(res_start, res_end + 1):
    vecs = []
    for label in dict_coil_per_residue:
        if res in dict_coil_per_residue[label]:
            vecs.append(dict_coil_per_residue[label][res])
    if vecs:
        concatenated_coil[res] = np.concatenate(vecs)

# -----------------------
# Concatenate CS
# -----------------------
concatenated_cs = {}
atom_types = ["N"]

for cs_file in CS_data:
    if os.path.exists(cs_file):
        df = read_cs_csv(cs_file)
        df_region, frame_cols = extract_region_cs(df, atom_types)

        for _, row in df_region.iterrows():
            res = row["resSeq"]
            cs_values = row[frame_cols].values.astype(float)

            if res not in concatenated_cs:
                concatenated_cs[res] = []

            concatenated_cs[res].append(cs_values)

for res in concatenated_cs:
    concatenated_cs[res] = np.concatenate(concatenated_cs[res])

# -----------------------
# Correlations
# -----------------------
correlations = {}
p_values = {}

residues_with_data = []
correlation_values = []
correlation_signed = []
p_value_list = []

for res in range(res_start, res_end + 1):
    if res in concatenated_coil and res in concatenated_cs:

        coil = concatenated_coil[res]
        cs = concatenated_cs[res]

        min_len = min(len(coil), len(cs))
        coil = coil[:min_len]
        cs = cs[:min_len]

        valid = ~np.isnan(cs)
        coil = coil[valid]
        cs = cs[valid]

        if len(coil) > 2 and np.std(cs) > 0:
            r, p = pearsonr(coil, cs)

            correlations[res] = r
            p_values[res] = p

            residues_with_data.append(res)
            correlation_values.append(abs(r))
            correlation_signed.append(r)
            p_value_list.append(p)

# -----------------------
# Scatter correlation plot (KEEP)
# -----------------------
if residues_with_data:

    cmap = mpl.colors.LinearSegmentedColormap.from_list(
        "white_orange_red", ["white", "orange", "red"]
    )

    fig, ax = plt.subplots(figsize=(15, 3))

    sc = ax.scatter(
        residues_with_data,
        correlation_values,
        c=correlation_values,
        cmap=cmap,
        s=120,
        edgecolors="black",
        linewidths=1.2,
        vmin=0,
        vmax=0.85
    )

    ax.set_xlabel("Residue")
    ax.set_ylabel("|Correlation|")
    ax.set_title("Coil vs N-CS correlation")

    tick_positions = [r for r in residues_with_data if r % 5 == 0]
    ax.set_xticks(tick_positions)

    ax.set_ylim(-0.05, 0.85)
    ax.set_xlim(res_start - 2, res_end)

    plt.colorbar(sc, ax=ax).set_label("|Pearson|")

    plt.tight_layout()
    plt.savefig(
        os.path.join(outdir, "pearson_correlation_coil_vs_N_CS_scatter.png"),
        dpi=300,
        bbox_inches="tight"
    )
    plt.close()

# -----------------------
# Top 5 individual residue plots (KEEP)
# -----------------------
print("Creating top5 plots...")

if len(residues_with_data) >= 5:

    top5_idx = np.argsort(correlation_values)[-5:][::-1]
    top5_res = [residues_with_data[i] for i in top5_idx]

    for res in top5_res:

        coil = concatenated_coil[res]
        cs = concatenated_cs[res]

        min_len = min(len(coil), len(cs))
        coil = coil[:min_len]
        cs = cs[:min_len]

        frames = np.arange(len(coil))

        fig, ax1 = plt.subplots(figsize=(12, 2))
        ax2 = ax1.twinx()

        ax1.scatter(frames, coil, color="grey", alpha=0.4, s=4)
        ax2.scatter(frames, cs, color="red", alpha=0.4, s=4)

        ax1.set_xlabel("Concatenated frames")
        ax1.set_ylabel("Coil (0/1)", color="grey")
        ax2.set_ylabel("N-CS (ppm)", color="red")

        ax1.set_ylim(-0.1, 1.1)
        ax1.set_yticks([0, 1])

        r = correlations[res]
        p = p_values[res]

        plt.title(
            f"Residue {res} (r={r:+.3f}, p={p:.4f})",
            weight="bold"
        )

        plt.tight_layout()
        plt.savefig(
            os.path.join(outdir, f"res{res}_coil_vs_N_CS.png"),
            dpi=300,
            bbox_inches="tight"
        )
        plt.close()

print("Done.")
