import os
import re
import glob
import pandas as pd
import matplotlib.pyplot as plt

"""
============================================================
 BEAST diffusion rate extractor (phylogenetic BD)
============================================================

Reads BEAST .log files under subdirectories and summarizes
the posterior of `coordinates.diffusionRate` after burn-in.

Supports two directory schemes via MODE:

1) Sampling Radius experiment
   - Subdirectories: r{radius}_sim{rep}   e.g., r10_sim3
   - Grouping key:  radius
   - X-axis label:  "Sampling Radius"
   - Figure name:   beast_diff_rate_per_radius.png

2) Maximal Dispersal Distance experiment
   - Subdirectories: dist{maxDist}_sim{rep}  e.g., dist10_sim3
   - Grouping key:   disp_dist_max
   - X-axis label:   "Max Dispersal Distance"
   - Figure name:    beast_diff_rate_per_maxDispDist.png
============================================================
"""

MODE = "maxdist"  # <-- set to "radius" or "maxdist"

if MODE == "radius":
    DIR_PATTERN = r"r(\d+)_sim(\d+)"
    KEY_NAME = "radius"
    XLABEL = "Sampling Radius"
    TITLE = "Estimated Diffusion Rate by Sampling Radius"
    FIG_NAME = "beast_diff_rate_per_radius.png"
elif MODE == "maxdist":
    DIR_PATTERN = r"dist(\d+)_sim(\d+)"
    KEY_NAME = "disp_dist_max"
    XLABEL = "Max Dispersal Distance"
    TITLE = "Estimated Diffusion Rate by Maximal Dispersal Distance"
    FIG_NAME = "beast_diff_rate_per_maxDispDist.png"
else:
    raise ValueError("MODE must be 'radius' or 'maxdist'")

def extract_and_save_dispersal_stats(
    burn_in_value=1_000_000,
    out_summary="beast_data.tsv",
    out_grouped="beast_data_grouped.tsv"
):
    results = []

    for directory in os.listdir("."):
        if not os.path.isdir(directory):
            continue

        m = re.match(DIR_PATTERN, directory)
        if not m:
            # not part of the experiment layout
            continue

        key_value, sim_num = map(int, m.groups())

        os.chdir(directory)
        try:
            log_files = glob.glob("*.log")
            if not log_files:
                continue

            # deterministic choice if multiple logs exist
            log_file = sorted(log_files, key=os.path.getmtime)[-1]

            df = pd.read_csv(log_file, comment="#", delimiter="\t")
            if "state" not in df.columns or "coordinates.diffusionRate" not in df.columns:
                # malformed log, skip
                continue

            # numeric coercion (in case of mixed types)
            df["state"] = pd.to_numeric(df["state"], errors="coerce")
            df["coordinates.diffusionRate"] = pd.to_numeric(
                df["coordinates.diffusionRate"], errors="coerce"
            )

            filtered_df = df[df["state"] > burn_in_value]

            # Delete .trees files to save space (best-effort)
            for tf in glob.glob("*.trees"):
                try:
                    os.remove(tf)
                except OSError:
                    pass

            if not filtered_df.empty:
                results.append({
                    KEY_NAME: key_value,
                    "simulation": sim_num,
                    "mean": filtered_df["coordinates.diffusionRate"].mean(),
                    "stdev": filtered_df["coordinates.diffusionRate"].std()
                })
        finally:
            os.chdir("..")

    global_df = pd.DataFrame(results)
    global_df.to_csv(out_summary, sep="\t", index=False)

    if global_df.empty:
        print("No valid BEAST results found; nothing to summarize/plot.")
        return pd.DataFrame()

    summary_df = (
        global_df.groupby(KEY_NAME)["mean"]
        .agg([
            ("mean_diffusion_rate", "mean"),
            ("median_diffusion_rate", "median"),
            ("quantile_2.5", lambda x: x.quantile(0.025)),
            ("quantile_97.5", lambda x: x.quantile(0.975))
        ])
        .reset_index()
    )
    summary_df.to_csv(out_grouped, sep="\t", index=False)
    return summary_df

summary_df = extract_and_save_dispersal_stats()

# Plot with quantile band if we have data
if not summary_df.empty:
    plt.figure(figsize=(10, 6))
    plt.plot(summary_df[KEY_NAME], summary_df["mean_diffusion_rate"],
             marker='o', label="Mean Estimated Diffusion Rate")

    plt.fill_between(
        summary_df[KEY_NAME],
        summary_df["quantile_2.5"],
        summary_df["quantile_97.5"],
        alpha=0.2,
        label="2.5% - 97.5% Quantile Range"
    )

    plt.xlabel(XLABEL)
    plt.ylabel("Estimated Diffusion Rate")
    plt.title(TITLE)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(FIG_NAME, dpi=300)
    # plt.show()
