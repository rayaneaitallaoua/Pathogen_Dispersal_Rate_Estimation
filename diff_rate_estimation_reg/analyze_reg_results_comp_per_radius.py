import glob
import os
import re
import pandas as pd

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
                parts = directory.split("_")
                radius = int(parts[0][1:])  # e.g., 'r10' -> 10
                sim_num = int(parts[1][3:])  # e.g., 'sim3' -> 3

                # Load iso file
                iso_file = glob.glob("*.ISO")
                if not iso_file:
                    continue

                try:
                    result = extract_slope_from_iso(iso_file[0])
                    results.append((radius, *result))

                except Exception as e:
                    print(f"Error reading {iso_file[0]} in {directory}: {e}")
                    continue

            finally:
                os.chdir("..")


    df = pd.DataFrame(results, columns=["radius", "slope", "ci_lower", "ci_upper"])

    df.to_csv("genepop_results.tsv", sep='\t', index=False)

    grouped = df.groupby("radius")["slope"].agg(["mean", "median", lambda x: x.quantile(0.025), lambda x: x.quantile(0.975)])
    grouped.columns = ["mean", "median", "quantile_2.5", "quantile_97.5"]
    grouped.to_csv("genepop_slope_summary.tsv", sep='\t')

aggregate_slopes()
