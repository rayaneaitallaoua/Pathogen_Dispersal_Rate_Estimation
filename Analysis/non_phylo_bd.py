import os
import re
from io import StringIO
import pandas as pd
import glob
import numpy as np

"""
============================================================
 Empirical diffusion rate extraction (non-phylogenetic BD)
============================================================

This script computes diffusion rates as var(X) / TMRCA
from GSpace and BEAST outputs. It supports two modes:

1) Sampling Radius experiment
   - Subdirectories: r{radius}_sim{rep}
   - Records grouped by: radius

2) Maximal Dispersal Distance experiment
   - Subdirectories: dist{maxDist}_sim{rep}
   - Records grouped by: max_distance

USAGE:
- Set MODE = "radius" or "maxdist" below.
- Output files:
    non_phylo_bd_all.tsv      (all replicate results)
    non_phylo_bd_grouped.tsv  (summary per radius or maxDist)
============================================================
"""

MODE = "maxdist"   # <-- change to "radius" or "maxdist"

# Define parsing patterns depending on the mode
if MODE == "radius":
    DIR_PATTERN = r"r(\d+)_sim(\d+)"
    KEY_NAME = "radius"
elif MODE == "maxdist":
    DIR_PATTERN = r"dist(\d+)_sim(\d+)"
    KEY_NAME = "max_distance"
else:
    raise ValueError("MODE must be either 'radius' or 'maxdist'")

def extract_empirical_diffusion_rates():
    records = []
    for directory in os.listdir("."):
        if not os.path.isdir(directory):
            continue
        match = re.match(DIR_PATTERN, directory)
        if not match:
            print(f"Skipping {directory}: name does not match {DIR_PATTERN}")
            continue

        key_value, sim_num = map(int, match.groups())
        os.chdir(directory)
        try:
            # Load GSpace summary
            gspace_file = glob.glob("*_seq_char.txt")
            if not gspace_file:
                os.chdir("..")
                continue
            gspace_df = pd.read_csv(gspace_file[0], sep=r'\s+', header=0)
            if "X" not in gspace_df.columns:
                print(f"Missing 'X' column in {gspace_file[0]} in {directory}")
                os.chdir("..")
                continue
            var_x = gspace_df["X"].var(ddof=0)  # population variance

            # Load BEAST log to get mean_age_root
            log_files = glob.glob("*.log")
            if not log_files:
                os.chdir("..")
                continue
            try:
                with open(log_files[0]) as f:
                    lines = [line for line in f if not line.startswith("Trace") and not line.startswith("#")]
                log_df = pd.read_csv(StringIO("".join(lines)), sep="\t", on_bad_lines="error")
            except Exception as e:
                print(f"Error reading {log_files[0]} in {directory}: {e}")
                os.chdir("..")
                continue

            if "age(root)" not in log_df.columns or "state" not in log_df.columns or log_df.empty:
                print(f"Invalid log file {log_files[0]} in {directory}")
                os.chdir("..")
                continue

            log_df["age(root)"] = pd.to_numeric(log_df["age(root)"], errors="coerce")
            burnin_index = int(0.1 * len(log_df))
            if burnin_index >= len(log_df):
                os.chdir("..")
                continue
            state_threshold = log_df["state"].iloc[burnin_index]
            log_df_post = log_df[log_df["state"] > state_threshold]

            if not log_df_post.empty:
                mean_age_root = log_df_post["age(root)"].mean()
                calc_diff_rate = (var_x / mean_age_root
                                  if pd.notna(mean_age_root) and mean_age_root > 0
                                  else np.nan)
                records.append({
                    KEY_NAME: key_value,
                    "simulation": sim_num,
                    "varX": var_x,
                    "mean_age_root": mean_age_root,
                    "calc_diff_rate": calc_diff_rate
                })
        finally:
            os.chdir("..")

    df = pd.DataFrame(records)
    df.to_csv("non_phylo_bd_all.tsv", sep="\t", index=False)

    if df.empty:
        print("No valid records found. Exiting.")
        return df

    grouped = (
        df.groupby(KEY_NAME)["calc_diff_rate"]
        .agg([
            ("mean", "mean"),
            ("median", "median"),
            ("quantile_2.5", lambda x: x.quantile(0.025)),
            ("quantile_97.5", lambda x: x.quantile(0.975)),
        ])
        .reset_index()
    )
    grouped.to_csv("non_phylo_bd_grouped.tsv", sep="\t", index=False)
    return df

if __name__ == "__main__":
    extract_empirical_diffusion_rates()
