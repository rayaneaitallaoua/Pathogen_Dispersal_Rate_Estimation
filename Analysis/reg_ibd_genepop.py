import glob
import os
import re
import pandas as pd

# using genepop IBD inverse slope

def extract_slope_from_iso(file_path):
    # Open the file and read all lines
    with open(file_path, 'r') as f:
        lines = f.readlines()
    # Iterate through lines to find the slope results block
    for i, line in enumerate(lines):
        if "ABC bootstrap results for SLOPE:" in line:
            # The slope result is expected on the next line
            target = lines[i + 2].strip()
            # Use regex to parse slope value and confidence intervals from the line
            match = re.search(r"([-+]?[0-9]*\.?[0-9]+) \[\s*([-+]?[0-9]*\.?[0-9]+)\s*,\s*([-+]?[0-9]*\.?[0-9]+)\s*\]", target)
            if match:
                slope, ci_low, ci_high = map(float, match.groups())
                # Return a tuple containing radius, slope, and confidence intervals
                return slope, ci_low, ci_high
    return None

def aggregate_slopes():

    results = []

    for directory in os.listdir("."):
        if os.path.isdir(directory):
            os.chdir(directory)
            try:
                # adapt pattern for radius
                match = re.match(r"dist(\d+)_sim(\d+)", directory)
                if not match:
                    print(f"Skipping directory {directory}: name does not match expected format.")
                    continue
                max_distance, sim_num = map(int, match.groups())

                # Load iso file
                iso_file = glob.glob("*.ISO")
                if not iso_file:
                    continue

                try:
                    result = extract_slope_from_iso(iso_file[0])
                    if result is None:
                        print(f"Skipping {directory}: could not parse slope from {iso_file[0]}.")
                        continue
                    slope, ci_low, ci_high = result

                    inv_slope = float("nan") if slope == 0 else 1.0 / slope
                    inv_ci_low = float("nan") if ci_low == 0 else 1.0 / ci_low
                    inv_ci_high = float("nan") if ci_high == 0 else 1.0 / ci_high

                    if not pd.isna(inv_ci_low) and not pd.isna(inv_ci_high) and inv_ci_low > inv_ci_high:
                        inv_ci_low, inv_ci_high = inv_ci_high, inv_ci_low
                    results.append((max_distance, slope, ci_low, ci_high, inv_slope, inv_ci_low, inv_ci_high))

                except Exception as e:
                    print(f"Error reading {iso_file[0]} in {directory}: {e}")
                    continue

            finally:
                os.chdir("..")

    df = pd.DataFrame(results, columns=["max_distance", "slope", "ci_lower", "ci_upper",
                                        "inv_slope", "inv_ci_lower", "inv_ci_upper"])

    df.to_csv("genepop_results.tsv", sep='\t', index=False)

    inv_results = []

    for _, row in df.iterrows():
        slope = row["slope"]
        ci_low = row["ci_lower"]
        ci_high = row["ci_upper"]
        max_distance = row["max_distance"]

        inv_slope = float("nan") if slope == 0 else 1.0 / slope
        inv_ci_low = float("nan") if ci_low == 0 else 1.0 / ci_low
        inv_ci_high = float("nan") if ci_high == 0 else 1.0 / ci_high

        if not pd.isna(inv_ci_low) and not pd.isna(inv_ci_high) and inv_ci_low > inv_ci_high:
            inv_ci_low, inv_ci_high = inv_ci_high, inv_ci_low

        inv_results.append((max_distance, inv_slope, inv_ci_low, inv_ci_high))

    df_inv = pd.DataFrame(inv_results, columns=["max_distance", "inv_slope", "inv_ci_lower", "inv_ci_upper"])
    df_inv.to_csv("genepop_inv_only.tsv", sep='\t', index=False)

    grouped = df.groupby("max_distance")["slope"].agg(["mean", "median", lambda x: x.quantile(0.025), lambda x: x.quantile(0.975)])
    grouped.columns = ["mean", "median", "quantile_2.5", "quantile_97.5"]
    grouped.to_csv("genepop_slope_summary.tsv", sep='\t')

    # Summary statistics for 1/slope
    grouped_inv = df.groupby("max_distance")["inv_slope"].agg(["mean", "median", lambda x: x.quantile(0.025), lambda x: x.quantile(0.975)])
    grouped_inv.columns = ["mean", "median", "quantile_2.5", "quantile_97.5"]
    grouped_inv.to_csv("genepop_invSlope_summary.tsv", sep='\t')


aggregate_slopes()
