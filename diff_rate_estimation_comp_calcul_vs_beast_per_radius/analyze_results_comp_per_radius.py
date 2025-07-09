import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
from io import StringIO

def extract_beast_diffusion_rates(burn_in_fraction=0.1):
    records = []
    for directory in os.listdir("."):
        if not os.path.isdir(directory):
            continue
        os.chdir(directory)
        try:
            parts = directory.split("_")
            radius = int(parts[0][1:])   # e.g., 'r10' -> 10
            sim_num = int(parts[1][3:])  # e.g., 'sim3' -> 3

            # Load BEAST log
            log_files = glob.glob("*.log")
            if not log_files:
                continue

            try:
                with open(log_files[0]) as f:
                    lines = [line for line in f if not line.startswith("Trace") and not line.startswith("#")]
                log_df = pd.read_csv(StringIO("".join(lines)), sep="\t",
                                     on_bad_lines="error")  # change to "warn" or "skip" to suppress
            except Exception as e:
                print(f"Error reading {log_files[0]} in {directory}: {e}")
                continue

            log_df["coordinates.diffusionRate"] = pd.to_numeric(log_df["coordinates.diffusionRate"], errors="coerce")
            burnin_index = int(burn_in_fraction * len(log_df))
            state_threshold = log_df["state"].iloc[burnin_index]
            log_df_post = log_df[log_df["state"] > state_threshold]

            # Delete .trees file to save space
            tree_files = glob.glob("*.trees")
            for tf in tree_files:
                os.remove(tf)

            if not log_df_post.empty:
                mean_diff_rate = log_df_post["coordinates.diffusionRate"].mean()
                records.append({
                    "radius": radius,
                    "simulation": sim_num,
                    "mean_diff_rate_beast": mean_diff_rate
                })
        finally:
            os.chdir("..")
    df = pd.DataFrame(records)
    df.to_csv("beast_diffusion_all.tsv", sep="\t", index=False)

    grouped = (
        df.groupby("radius")["mean_diff_rate_beast"]
        .agg([
            ("mean", "mean"),
            ("median", "median"),
            ("quantile_2.5", lambda x: x.quantile(0.025)),
            ("quantile_97.5", lambda x: x.quantile(0.975)),
        ])
        .reset_index()
    )
    grouped.to_csv("beast_diffusion_grouped.tsv", sep="\t", index=False)
    return df

def extract_empirical_diffusion_rates():
    records = []
    for directory in os.listdir("."):
        if not os.path.isdir(directory):
            continue
        os.chdir(directory)
        try:
            parts = directory.split("_")
            radius = int(parts[0][1:])   # e.g., 'r10' -> 10
            sim_num = int(parts[1][3:])  # e.g., 'sim3' -> 3

            # Load GSpace summary
            gspace_file = glob.glob("*_seq_char.txt")
            if not gspace_file:
                continue
            gspace_df = pd.read_csv(gspace_file[0], sep=r'\s+', header=0)
            var_x = gspace_df["X"].var()

            # Load BEAST log to get mean_age_root
            log_files = glob.glob("*.log")
            if not log_files:
                continue

            try:
                with open(log_files[0]) as f:
                    lines = [line for line in f if not line.startswith("Trace") and not line.startswith("#")]
                log_df = pd.read_csv(StringIO("".join(lines)), sep="\t",
                                     on_bad_lines="error")  # change to "warn" or "skip" to suppress
            except Exception as e:
                print(f"Error reading {log_files[0]} in {directory}: {e}")
                continue
            log_df["age(root)"] = pd.to_numeric(log_df["age(root)"], errors="coerce")
            burnin_index = int(0.1 * len(log_df))
            state_threshold = log_df["state"].iloc[burnin_index]
            log_df_post = log_df[log_df["state"] > state_threshold]
            if not log_df_post.empty:
                mean_age_root = log_df_post["age(root)"].mean()
                calc_diff_rate = var_x / mean_age_root if mean_age_root > 0 else None
                records.append({
                    "radius": radius,
                    "simulation": sim_num,
                    "varX": var_x,
                    "mean_age_root": mean_age_root,
                    "calc_diff_rate": calc_diff_rate
                })
        finally:
            os.chdir("..")
    df = pd.DataFrame(records)
    df.to_csv("empirical_diffusion_all.tsv", sep="\t", index=False)

    grouped = (
        df.groupby("radius")["calc_diff_rate"]
        .agg([
            ("mean", "mean"),
            ("median", "median"),
            ("quantile_2.5", lambda x: x.quantile(0.025)),
            ("quantile_97.5", lambda x: x.quantile(0.975)),
        ])
        .reset_index()
    )
    grouped.to_csv("empirical_diffusion_grouped.tsv", sep="\t", index=False)
    return df

def plot_diffusion_rate_comparison(beast_df, empirical_df):
    beast_grouped = pd.read_csv("beast_diffusion_grouped.tsv", sep="\t")
    empirical_grouped = pd.read_csv("empirical_diffusion_grouped.tsv", sep="\t")

    plt.figure(figsize=(10, 6))
    plt.plot(empirical_grouped["radius"], empirical_grouped["mean"], marker='o', label="Var(X)/T Estimate")
    plt.plot(beast_grouped["radius"], beast_grouped["mean"], marker='s', label="BEAST Estimate")
    plt.fill_between(
        empirical_grouped["radius"],
        empirical_grouped["quantile_2.5"],
        empirical_grouped["quantile_97.5"],
        color='blue',
        alpha=0.2,
        label="Var(X)/T 95% CI"
    )
    plt.fill_between(
        beast_grouped["radius"],
        beast_grouped["quantile_2.5"],
        beast_grouped["quantile_97.5"],
        color='orange',
        alpha=0.2,
        label="BEAST 95% CI"
    )
    plt.xlabel("Radius")
    plt.ylabel("Diffusion Rate")
    plt.title("Comparison of Diffusion Rate Estimates")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("diffusion_rate_comparison.png", dpi=300)

beast_df = extract_beast_diffusion_rates()
empirical_df = extract_empirical_diffusion_rates()
#plot_diffusion_rate_comparison(beast_df, empirical_df)