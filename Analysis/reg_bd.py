import os
import re
import glob
import math
import itertools
import pandas as pd
from Bio import SeqIO
from scipy.stats import linregress
from scipy.spatial.distance import euclidean

"""
============================================================
 Reverse regression (BD): geo_d^2 ~ gen_d
============================================================

Computes pairwise Jukes–Cantor genetic distances (gen_d) and
squared Euclidean geographic distances (geo_d^2) from FASTA
headers, then fits the reverse regression:

    linregress(gen_d, geo_d2)

Supports two directory schemes and header coordinate layouts:

1) Sampling Radius experiment
   - Subdirectories: r{radius}_sim{rep}      (e.g., r10_sim3)
   - FASTA header coords at indices [9], [10]
   - Grouping key:  "radius"

2) Maximal Dispersal Distance experiment
   - Subdirectories: dist{maxDist}_sim{rep}  (e.g., dist10_sim3)
   - FASTA header coords at indices [10], [11]
   - Grouping key:  "disp_dist_max"

USAGE:
- Set MODE = "radius" or "maxdist" below; the script adapts
  dir parsing, coordinate extraction, grouping, and outputs.

Outputs:
- regression_bd_slopes.tsv   (per replicate)
- regression_bd_summary.tsv  (summary per radius / maxDist)
============================================================
"""

MODE = "maxdist"  # <-- set to "radius" or "maxdist"

if MODE == "radius":
    DIR_PATTERN = r"r(\d+)_sim(\d+)"
    KEY_NAME = "radius"
    COORD_IDX = (9, 10)     # x at [9], y at [10]
elif MODE == "maxdist":
    DIR_PATTERN = r"dist(\d+)_sim(\d+)"
    KEY_NAME = "disp_dist_max"
    COORD_IDX = (10, 11)    # x at [10], y at [11]
else:
    raise ValueError("MODE must be 'radius' or 'maxdist'")

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
        if len(header_parts) <= max(COORD_IDX):
            print(f"Invalid header format in {rec.id}")
            continue
        try:
            x = float(header_parts[COORD_IDX[0]])
            y = float(header_parts[COORD_IDX[1]])
        except ValueError:
            print(f"Invalid coordinate format in {rec.id}")
            continue
        coords[rec.id] = (x, y)
        sequences[rec.id] = str(rec.seq)
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
        return (float("nan"),) * 4

    df = pd.DataFrame(pairs, columns=["geo_d2", "gen_d"])
    if df["geo_d2"].nunique() < 2 or df["gen_d"].nunique() < 2:
        print("Warning: insufficient variation — skipping regression.")
        return (float("nan"),) * 4

    slope, _, r_value, p_value, std_err = linregress(df["gen_d"], df["geo_d2"])
    return slope, r_value**2, p_value, std_err

results = []

for directory in os.listdir("."):
    if not os.path.isdir(directory):
        continue

    m = re.match(DIR_PATTERN, directory)
    if not m:
        continue

    key_value, sim_num = map(int, m.groups())

    os.chdir(directory)
    try:
        fa_files = glob.glob("*.fa")
        if not fa_files:
            print(f"No .fa file in {directory}")
            continue

        # pick a most recent FASTA if multiple exist
        fasta_path = sorted(fa_files, key=os.path.getmtime)[-1]
        sequences, coords = load_sequences_and_coords(fasta_path)

        if len(sequences) < 2:
            print(f"Too few valid sequences in {fasta_path}")
            continue

        slope, r2, pval, stderr = compute_pairwise_regression(sequences, coords)
        results.append((key_value, sim_num, slope, r2, pval, stderr))
    finally:
        os.chdir("..")

df_out = pd.DataFrame(results, columns=[
    KEY_NAME, "sim", "slope", "r2", "pval", "stderr"
])
df_out.to_csv("regression_bd_slopes.tsv", sep="\t", index=False)

summary_df = df_out.groupby(KEY_NAME).agg({
    "slope": ["mean", "median",
              lambda x: x.quantile(0.025),
              lambda x: x.quantile(0.975)]
})

summary_df.columns = ["slope_mean", "slope_median", "slope_q2.5", "slope_q97.5"]
summary_df = summary_df.reset_index()
summary_df.to_csv("regression_bd_summary.tsv", sep="\t", index=False)
