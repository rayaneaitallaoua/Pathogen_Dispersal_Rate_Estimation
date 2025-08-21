import os
import re
import glob
import pandas as pd
import math

"""
============================================================
 Genepop IBD results aggregator (slope and 1/slope)
============================================================

Parses *.ISO outputs to extract the a(r) slope and its
bootstrap CI, computes 1/slope and its inverted CI, and
summarizes results by replicate and by group.

Supports two directory schemes via MODE:

1) Sampling Radius experiment
   - Subdirectories: r{radius}_sim{rep}   (e.g., r10_sim3)
   - Grouping key:  "radius"

2) Maximal Dispersal Distance experiment
   - Subdirectories: dist{maxDist}_sim{rep}  (e.g., dist10_sim3)
   - Grouping key:   "max_distance"

USAGE:
- Set MODE = "radius" or "maxdist" below; the script adapts
  directory matching and the grouping key automatically.

Outputs:
- genepop_results.tsv          (per replicate: slope & 1/slope)
- genepop_inv_only.tsv         (per replicate: 1/slope only)
- genepop_slope_summary.tsv    (summary by group for slope)
- genepop_invSlope_summary.tsv (summary by group for 1/slope)
============================================================
"""

MODE = "maxdist"  # <-- set to "radius" or "maxdist"

if MODE == "radius":
    DIR_PATTERN = r"r(\d+)_sim(\d+)"
    KEY_NAME = "radius"
elif MODE == "maxdist":
    DIR_PATTERN = r"dist(\d+)_sim(\d+)"
    KEY_NAME = "max_distance"
else:
    raise ValueError("MODE must be 'radius' or 'maxdist'")

# robust numeric pattern (handles scientific notation)
_NUM = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"
_SLOPE_CI_RE = re.compile(rf"({_NUM})\s*\[\s*({_NUM})\s*,\s*({_NUM})\s*\]")

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
        if not os.path.isdir(directory):
            continue

        m = re.match(DIR_PATTERN, directory)
        if not m:
            # Not a target directory for the chosen MODE
            continue

        key_value, sim_num = map(int, m.groups())

        os.chdir(directory)
        try:
            iso_files = glob.glob("*.ISO")
            if not iso_files:
                continue

            # pick most recent ISO if multiple
            iso_file = sorted(iso_files, key=os.path.getmtime)[-1]

            try:
                parsed = extract_slope_from_iso(iso_file)
            except Exception as e:
                print(f"Error parsing {iso_file} in {directory}: {e}")
                continue

            if parsed is None:
                print(f"Skipping {directory}: could not parse slope from {iso_file}.")
                continue

            slope, ci_low, ci_high = parsed

            inv_slope = math.nan if slope == 0 else 1.0 / slope
            inv_ci_low = math.nan if ci_low == 0 else 1.0 / ci_low
            inv_ci_high = math.nan if ci_high == 0 else 1.0 / ci_high

            # After inversion, bounds may flip; reorder if needed
            if not pd.isna(inv_ci_low) and not pd.isna(inv_ci_high) and inv_ci_low > inv_ci_high:
                inv_ci_low, inv_ci_high = inv_ci_high, inv_ci_low

            results.append(
                (key_value, sim_num, slope, ci_low, ci_high, inv_slope, inv_ci_low, inv_ci_high)
            )

        finally:
            os.chdir("..")

    # Per-replicate table
    df = pd.DataFrame(
        results,
        columns=[KEY_NAME, "simulation", "slope", "ci_lower", "ci_upper",
                 "inv_slope", "inv_ci_lower", "inv_ci_upper"]
    )
    df.to_csv("genepop_results.tsv", sep="\t", index=False)

    # 1/slope only per replicate (optional convenience)
    df_inv = df[[KEY_NAME, "inv_slope", "inv_ci_lower", "inv_ci_upper"]].copy()
    df_inv.to_csv("genepop_inv_only.tsv", sep="\t", index=False)

    if df.empty:
        print("No valid ISO results found; summaries not generated.")
        return

    # Summary statistics for slope
    grouped = (
        df.groupby(KEY_NAME)["slope"]
          .agg([
              ("mean", "mean"),
              ("median", "median"),
              ("quantile_2.5", lambda x: x.quantile(0.025)),
              ("quantile_97.5", lambda x: x.quantile(0.975)),
          ])
          .reset_index()
    )
    grouped.to_csv("genepop_slope_summary.tsv", sep="\t", index=False)

    # Summary statistics for 1/slope
    grouped_inv = (
        df.groupby(KEY_NAME)["inv_slope"]
          .agg([
              ("mean", "mean"),
              ("median", "median"),
              ("quantile_2.5", lambda x: x.quantile(0.025)),
              ("quantile_97.5", lambda x: x.quantile(0.975)),
          ])
          .reset_index()
    )
    grouped_inv.to_csv("genepop_invSlope_summary.tsv", sep="\t", index=False)

if __name__ == "__main__":
    aggregate_slopes()
