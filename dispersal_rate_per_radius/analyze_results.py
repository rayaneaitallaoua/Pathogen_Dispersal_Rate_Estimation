import os
import glob
import pandas as pd
import matplotlib.pyplot as plt

burn_in_value = 1000000
results = []

def extract_dispersal_rate(burn_in):
    log_files = glob.glob("*.log")
    if not log_files:
        return None, None
    log_file_path = log_files[0]

    df = pd.read_csv(log_file_path, comment='#', delimiter='\t')
    filtered_df = df[df["state"] > burn_in]

    return filtered_df["coordinates.diffusionRate"].mean(), filtered_df["coordinates.diffusionRate"].std()

for directory in os.listdir("."):
    if not os.path.isdir(directory):
        continue

    os.chdir(directory)

    try:
        parts = directory.split("_")
        radius = int(parts[0][1:])   # e.g., 'r10' -> 10
        sim_num = int(parts[1][3:])  # e.g., 'sim3' -> 3

        mean, stdev = extract_dispersal_rate(burn_in_value)
        if mean is not None:
            results.append({
                "radius": radius,
                "simulation": sim_num,
                "mean": mean,
                "stdev": stdev
            })
    finally:
        os.chdir("..")

# Convert list of results to DataFrame
global_df = pd.DataFrame(results)
global_df.to_csv("dispersal_summary.tsv", sep="\t", index=False)

# Group by radius and compute summary statistics
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

# Save as tab-separated file
summary_df.to_csv("dispersal_summary_grouped.tsv", sep="\t", index=False)

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
