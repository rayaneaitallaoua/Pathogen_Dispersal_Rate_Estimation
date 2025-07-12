import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy.stats import f_oneway

# Load data
beast_df = pd.read_csv("all_beast_data.tsv", sep="\t")
emp_df = pd.read_csv("all_data_emp.tsv", sep="\t")
reg_df = pd.read_csv("all_data_reg.tsv", sep="\t")

# Standardize columns
beast_df = beast_df.rename(columns={"mean_diff_rate_beast": "estimation"})
emp_df = emp_df.rename(columns={"calc_diff_rate": "estimation"})
reg_df = reg_df.rename(columns={"slope": "estimation"})

beast_df = beast_df[["radius", "sample_size", "method", "simulation", "estimation"]]
emp_df = emp_df[["radius", "sample_size", "method", "simulation", "estimation"]]

reg_df["simulation"] = reg_df.groupby(["radius", "sample_size"]).cumcount() + 1
reg_df = reg_df[["radius", "sample_size", "method", "simulation", "estimation"]]

# Combine all
combined_df = pd.concat([beast_df, emp_df, reg_df], ignore_index=True)
combined_df = combined_df[["radius", "simulation", "sample_size", "method", "estimation"]]
combined_df = combined_df.sort_values(by=["method", "radius", "simulation"], ignore_index=True)

# Convert radius and sample_size to categorical for ANOVA
combined_df["radius"] = combined_df["radius"].astype("category")
combined_df["sample_size"] = combined_df["sample_size"].astype("category")

output = []

# 1-Way ANOVA across all methods
output.append("=== 1-Way ANOVA: estimation ~ method (all data) ===")
output.append("# Question: Is there a difference between estimation methods across all radii and sample sizes?")
method_groups = [g["estimation"].values for _, g in combined_df.groupby("method")]
stat, p = f_oneway(*method_groups)
output.append(f"F = {stat:.4f}, p = {p:.4e}")
output.append("# ✅ Significant difference between methods" if p < 0.05 else "# ❌ No significant difference between methods")
output.append("")

# 1-Way ANOVA: estimation ~ radius (per method)
output.append("=== 1-Way ANOVA: estimation ~ radius (within each method) ===")
output.append("# Question: Is the variation of estimation between radii significant for each method?")
for method in combined_df["method"].unique():
    subset = combined_df[combined_df["method"] == method]
    groups = [g["estimation"].values for _, g in subset.groupby("radius")]
    if all(len(g) > 1 for g in groups):
        stat, p = f_oneway(*groups)
        output.append(f"Method: {method} --> F = {stat:.4f}, p = {p:.4e}")
        output.append("# ✅ Significant variation between radii" if p < 0.05 else "# ❌ No significant variation between radii")
output.append("")

# 1-Way ANOVA: estimation ~ sample_size (per method)
output.append("=== 1-Way ANOVA: estimation ~ sample_size (within each method) ===")
output.append("# Question: Is the variation of estimation between sample sizes significant for each method?")
for method in combined_df["method"].unique():
    subset = combined_df[combined_df["method"] == method]
    groups = [g["estimation"].values for _, g in subset.groupby("sample_size")]
    if all(len(g) > 1 for g in groups):
        stat, p = f_oneway(*groups)
        output.append(f"Method: {method} --> F = {stat:.4f}, p = {p:.4e}")
        output.append("# ✅ Significant variation between sample sizes" if p < 0.05 else "# ❌ No significant variation between sample sizes")
output.append("")

# Save output
with open("stats.txt", "w") as f:
    f.writelines(line + "\n" for line in output)