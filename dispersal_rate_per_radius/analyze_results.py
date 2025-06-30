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
            radius = int(parts[0][1:])   # e.g., 'r10' -> 10
            sim_num = int(parts[1][3:])  # e.g., 'sim3' -> 3

            log_files = glob.glob("*.log")
            if not log_files:
                continue

            df = pd.read_csv(log_files[0], comment='#', delimiter='\t')
            filtered_df = df[df["state"] > burn_in_value]

            if not filtered_df.empty:
                results.append({
                    "radius": radius,
                    "simulation": sim_num,
                    "mean": filtered_df["coordinates.diffusionRate"].mean(),
                    "stdev": filtered_df["coordinates.diffusionRate"].std()
                })
        finally:
            os.chdir("..")

    global_df = pd.DataFrame(results)
    global_df.to_csv(out_summary, sep="\t", index=False)

    summary_df = (
        global_df.groupby("radius")["mean"]
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
plt.plot(summary_df["radius"], summary_df["mean_dispersal_rate"], marker='o', label="Mean Dispersal Rate")

# Add shaded area for quantile range
plt.fill_between(
    summary_df["radius"],
    summary_df["quantile_2.5"],
    summary_df["quantile_97.5"],
    color='blue',
    alpha=0.2,
    label="2.5% - 97.5% Quantile Range"
)

plt.xlabel("Radius")
plt.ylabel("Dispersal Rate")
plt.title("Dispersal Rate by Radius with Quantile Range")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

def extract_empirical_dispersal():
    empirical_dispersal_data = []

    for root, dirs, files in os.walk("."):
        for file in files:
            if file.endswith("_GSpace_param_summary.txt"):
                full_path = os.path.join(root, file)
                parts = os.path.normpath(full_path).split(os.sep)
                try:
                    radius = int([p for p in parts if p.startswith("r")][0][1:])
                    sim_num = int([p for p in parts if p.startswith("sim")][0][3:])
                except (IndexError, ValueError):
                    continue

                with open(full_path, "r") as f:
                    lines = f.readlines()
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
                                "radius": radius,
                                "simulation": sim_num,
                                "step": int(step),
                                "freqX": float(freqX),
                                "freqY": float(freqY)
                            })

    empirical_df = pd.DataFrame(empirical_dispersal_data)
    empirical_df.to_csv("empirical_dispersal_all.tsv", sep="\t", index=False)

    empirical_summary = (
        empirical_df.groupby(["radius", "step"])[["freqX", "freqY"]]
        .mean()
        .reset_index()
    )
    empirical_summary.to_csv("empirical_dispersal_grouped.tsv", sep="\t", index=False)

extract_empirical_dispersal()