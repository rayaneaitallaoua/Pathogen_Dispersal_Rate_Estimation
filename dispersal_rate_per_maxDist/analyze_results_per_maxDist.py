import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

def extract_and_save_dispersal_stats(burn_in_value=1000000,
                                     out_summary="dispersal_summary.tsv",
                                     out_grouped="dispersal_summary_grouped.tsv"):
    results = []

    for directory in os.listdir("."):
        if not os.path.isdir(directory):
            continue

        os.chdir(directory)
        try:
            parts = directory.split("_")
            disp_dist_max = int(parts[0][4:])   # e.g., 'dist10' -> 10
            sim_num = int(parts[1][3:])  # e.g., 'sim3' -> 3

            log_files = glob.glob("*.log")
            if not log_files:
                continue

            df = pd.read_csv(log_files[0], comment='#', delimiter='\t')
            filtered_df = df[df["state"] > burn_in_value]

            # Delete .trees file to save space
            tree_files = glob.glob("*.trees")
            for tf in tree_files:
                os.remove(tf)

            if not filtered_df.empty:
                results.append({
                    "disp_dist_max": disp_dist_max,
                    "simulation": sim_num,
                    "mean": filtered_df["coordinates.diffusionRate"].mean(),
                    "stdev": filtered_df["coordinates.diffusionRate"].std()
                })
        finally:
            os.chdir("..")

    global_df = pd.DataFrame(results)
    global_df.to_csv(out_summary, sep="\t", index=False)

    summary_df = (
        global_df.groupby("disp_dist_max")["mean"]
        .agg([
            ("mean_dispersal_rate", "mean"),
            ("median_dispersal_rate", "median"),
            ("quantile_2.5", lambda x: x.quantile(0.025)),
            ("quantile_97.5", lambda x: x.quantile(0.975))
        ])
        .reset_index()
    )
    summary_df.to_csv(out_grouped, sep="\t", index=False)
    return summary_df

summary_df = extract_and_save_dispersal_stats()

# Plot with error bars (quantile range)
plt.figure(figsize=(10, 6))
plt.plot(summary_df["disp_dist_max"], summary_df["mean_dispersal_rate"], marker='o', label="Mean Estimated Dispersal Rate")

# Add shaded area for quantile range
plt.fill_between(
    summary_df["disp_dist_max"],
    summary_df["quantile_2.5"],
    summary_df["quantile_97.5"],
    color='blue',
    alpha=0.2,
    label="2.5% - 97.5% Quantile Range"
)

plt.xlabel("Max Dispersal Distance")
plt.ylabel("Estimated Dispersal Rate")
plt.title("Estimated Dispersal Rate by Maximal Dispersal Distance")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("beast_dispersal_rate_per_dist.png",dpi=300)

def extract_empirical_dispersal(full_output="empirical_dispersal_all.tsv",
                                grouped_output="empirical_dispersal_grouped.tsv",
                                moment_output="empirical_dispersal_moments_all.tsv",
                                moment_grouped_output="empirical_dispersal_moments_grouped.tsv"):

    empirical_dispersal_data = []
    empirical_moments_data = []

    for directory in os.listdir("."):
        if not os.path.isdir(directory):
            continue

        os.chdir(directory)
        try:
            parts = directory.split("_")
            disp_dist_max = int(parts[0][4:])   # e.g., 'dist10' -> 10
            sim_num = int(parts[1][3:])  # e.g., 'sim3' -> 3

            for file in os.listdir("."):
                if file.endswith("_GSpace_param_summary.txt"):
                    with open(file, "r") as f:
                        lines = f.readlines()

                    # Extract dispersal distribution
                    start_idx = None
                    for i, line in enumerate(lines):
                        if line.strip().startswith("Empirical dispersal distribution"):
                            start_idx = i + 1
                            break

                    if start_idx is not None:
                        for j in range(start_idx, len(lines)):
                            line = lines[j].strip()
                            if not line or line.startswith("#") or line.startswith("%%"):
                                break
                            parts = line.split()
                            if len(parts) != 3:
                                continue
                            step, freqX, freqY = parts
                            empirical_dispersal_data.append({
                                "disp_dist_max": disp_dist_max,
                                "simulation": sim_num,
                                "step": int(step),
                                "freqX": float(freqX),
                                "freqY": float(freqY)
                            })

                    # Extract dispersal moments
                    meanX = meanY = varX = varY = skewX = skewY = kurtX = kurtY = None
                    for line in lines:
                        if line.startswith("Empirical dispersal mean"):
                            meanX, meanY = map(float, line.split(":")[1].strip().split())
                        elif line.startswith("Empirical dispersal sigma2"):
                            varX, varY = map(float, line.split(":")[1].strip().split())
                        elif line.startswith("Empirical dispersal skewness"):
                            skewX, skewY = map(float, line.split(":")[1].strip().split())
                        elif line.startswith("Empirical dispersal kurtosis"):
                            kurtX, kurtY = map(float, line.split(":")[1].strip().split())

                    empirical_moments_data.append({
                        "disp_dist_max": disp_dist_max,
                        "simulation": sim_num,
                        "meanX": meanX,
                        "meanY": meanY,
                        "varX": varX,
                        "varY": varY,
                        "skewX": skewX,
                        "skewY": skewY,
                        "kurtX": kurtX,
                        "kurtY": kurtY
                    })
        finally:
            os.chdir("..")

    empirical_df = pd.DataFrame(empirical_dispersal_data)
    empirical_df.to_csv(full_output, sep="\t", index=False)

    empirical_summary = (
        empirical_df.groupby(["disp_dist_max", "step"])[["freqX", "freqY"]]
        .mean()
        .reset_index()
    )
    empirical_summary.to_csv(grouped_output, sep="\t", index=False)

    moments_df = pd.DataFrame(empirical_moments_data)
    moments_df.to_csv(moment_output, sep="\t", index=False)

    grouped_moments = (
        moments_df.groupby("disp_dist_max")[["meanX", "meanY", "varX", "varY", "skewX", "skewY", "kurtX", "kurtY"]]
        .mean()
        .reset_index()
    )
    grouped_moments.to_csv(moment_grouped_output, sep="\t", index=False)

    # Plot each moment parameter (mean, var, skew, kurt) in a separate subplot with X and Y
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    axs = axs.flatten()
    parameters = [("meanX", "meanY", "Mean"), ("varX", "varY", "Variance"),
                  ("skewX", "skewY", "Skewness"), ("kurtX", "kurtY", "Kurtosis")]

    for i, (x_param, y_param, title) in enumerate(parameters):
        axs[i].plot(grouped_moments["disp_dist_max"], grouped_moments[x_param], label="X", color="blue")
        axs[i].plot(grouped_moments["disp_dist_max"], grouped_moments[y_param], label="Y", color="orange")
        axs[i].set_title(f"{title} by maximum dispersal distance")
        axs[i].set_xlabel("Max Dispersal Distance")
        axs[i].set_ylabel(title)
        axs[i].legend()
        axs[i].grid(True)

    plt.tight_layout()
    plt.savefig("empirical_dispersal_moments_by_dist.png",dpi=300)

extract_empirical_dispersal()