import os
import glob
import math
import itertools
import pandas as pd
from Bio import SeqIO
from scipy.stats import linregress
from scipy.spatial.distance import euclidean

def jukes_cantor_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be the same length")
    differences = sum(c1 != c2 for c1, c2 in zip(seq1, seq2))
    p = differences / len(seq1)
    if p >= 0.75:
        return float("nan")
    try:
        return -3/4 * math.log(1 - (4/3)*p)
    except ValueError:
        return float("nan")

def load_sequences_and_coords(fasta_path):
    sequences = {}
    coords = {}
    records = list(SeqIO.parse(fasta_path, "fasta"))
    for rec in records[1:]:  # skip ancestral
        header_parts = rec.id.split("_")
        if len(header_parts) < 12:
            print(f"Invalid header format in {rec.id}")
            continue
        try:
            # for radius
#            x = float(header_parts[9])
#            y = float(header_parts[10])

            # for max disp dist
            x = float(header_parts[10])
            y = float(header_parts[11])
            coords[rec.id] = (x, y)
            sequences[rec.id] = str(rec.seq)
        except ValueError:
            print(f"Invalid coordinate format in {rec.id}")
            continue
    return sequences, coords

def compute_pairwise_regression(sequences, coords):
    pairs = []
    for id1, id2 in itertools.combinations(sequences.keys(), 2):
        if id1 not in coords or id2 not in coords:
            continue
        geo_d2 = euclidean(coords[id1], coords[id2]) ** 2
        gen_d = jukes_cantor_distance(sequences[id1], sequences[id2])
        if not math.isnan(gen_d):
            pairs.append((geo_d2, gen_d))
    if len(pairs) < 2:
        return (float("nan"),)*4

    df = pd.DataFrame(pairs, columns=["geo_d2", "gen_d"])
    if df["geo_d2"].nunique() < 2 or df["gen_d"].nunique() < 2:
        print("Warning: insufficient variation — skipping regression.")
        return (float("nan"),)*4

    # regression (x, y)
    slope, _, r_value, p_value, std_err = linregress(df["gen_d"], df["geo_d2"])

    return slope, r_value**2, p_value, std_err

results = []

for directory in os.listdir("."):
    if not os.path.isdir(directory):
        continue

    os.chdir(directory)
    try:
        parts = directory.split("_")
        try:
#            radius = int(parts[0][1:])   # e.g., 'r10' -> 10
            disp_dist_max = int(parts[0][4:])   # e.g., 'dist10' -> 10
            sim_num = int(parts[1][3:])         # e.g., 'sim3' -> 3
        except (IndexError, ValueError):
            print(f"Skipping directory with unexpected name: {directory}")
            os.chdir("..")
            continue

        fa_files = glob.glob("*.fa")
        if not fa_files:
            print(f"No .fa file in {directory}")
            os.chdir("..")
            continue

        fasta_path = fa_files[0]
        sequences, coords = load_sequences_and_coords(fasta_path)

        if len(sequences) < 2:
            print(f"Too few valid sequences in {fasta_path}")
            os.chdir("..")
            continue

        slope, r2, pval, stderr = compute_pairwise_regression(sequences, coords)

        results.append((disp_dist_max, sim_num, slope, r2, pval, stderr))
    finally:
        os.chdir("..")

df_out = pd.DataFrame(results, columns=[
    "disp_dist_max", "sim",
    "slope", "r2", "pval", "stderr"
])

df_out.to_csv("regression_bd_slopes.tsv", sep="\t", index=False)

# Group by Maximal Dispersal Distance and compute summary stats
summary_df = df_out.groupby("disp_dist_max").agg({
    "slope": ["mean", "median", lambda x: x.quantile(0.025), lambda x: x.quantile(0.975)]
})

# Rename columns
summary_df.columns = [
    "slope_mean", "slope_median", "slope_q2.5", "slope_q97.5"
]
summary_df = summary_df.reset_index()

# Save summary
summary_df.to_csv("regression_bd_summary.tsv", sep="\t", index=False)
