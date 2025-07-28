import os
import re
from io import StringIO
import pandas as pd
import glob

# using var(x) / TMRCA
# with x being positions from GSpace simulation
# and TMRCA being root.age() mean from BEAST

def extract_empirical_diffusion_rates():
    records = []
    for directory in os.listdir("."):
        if not os.path.isdir(directory):
            continue
        os.chdir(directory)
        try:
            match = re.match(r"dist(\d+)_sim(\d+)", directory)
            if not match:
                print(f"Skipping directory {directory}: name does not match expected format.")
                os.chdir("..")
                continue
            max_distance, sim_num = map(int, match.groups())

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
            var_x = gspace_df["X"].var()

            # Load BEAST log to get mean_age_root
            log_files = glob.glob("*.log")
            if not log_files:
                os.chdir("..")
                continue

            try:
                with open(log_files[0]) as f:
                    lines = [line for line in f if not line.startswith("Trace") and not line.startswith("#")]
                log_df = pd.read_csv(StringIO("".join(lines)), sep="\t",
                                     on_bad_lines="error")  # change to "warn" or "skip" to suppress
            except Exception as e:
                print(f"Error reading {log_files[0]} in {directory}: {e}")
                os.chdir("..")
                continue
            if "age(root)" not in log_df.columns or "state" not in log_df.columns:
                print(f"Missing 'age(root)' or 'state' columns in {log_files[0]} in {directory}")
                os.chdir("..")
                continue
            log_df["age(root)"] = pd.to_numeric(log_df["age(root)"], errors="coerce")
            burnin_index = int(0.1 * len(log_df))
            state_threshold = log_df["state"].iloc[burnin_index]
            log_df_post = log_df[log_df["state"] > state_threshold]
            if not log_df_post.empty:
                mean_age_root = log_df_post["age(root)"].mean()
                calc_diff_rate = var_x / mean_age_root if mean_age_root > 0 else None
                records.append({
                    "max_distance": max_distance,
                    "simulation": sim_num,
                    "varX": var_x,
                    "mean_age_root": mean_age_root,
                    "calc_diff_rate": calc_diff_rate
                })
        finally:
            os.chdir("..")
    df = pd.DataFrame(records)
    df.to_csv("non_phylo_bd_all.tsv", sep="\t", index=False)

    grouped = (
        df.groupby("max_distance")["calc_diff_rate"]
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

extract_empirical_diffusion_rates()