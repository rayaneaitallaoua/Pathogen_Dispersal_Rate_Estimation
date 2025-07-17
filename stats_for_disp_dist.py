import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy.stats import f_oneway

# Load data
beast_data = pd.read_csv("/Users/ayoubrayaneaitallaoua/Desktop/IFB_analysis_result_archive/various_max_disp_dist/500_500_lattice/beast_data.tsv", sep="\t")
bd_data = pd.read_csv("/Users/ayoubrayaneaitallaoua/Desktop/IFB_analysis_result_archive/various_max_disp_dist/500_500_lattice/empirical_diffusion_all.tsv", sep="\t")
genepop_data = pd.read_csv("/Users/ayoubrayaneaitallaoua/Desktop/IFB_analysis_result_archive/various_max_disp_dist/500_500_lattice/genepop_results.tsv", sep="\t")

# Clean and rename
beast_data = beast_data.rename(columns={"disp_dist_max": "max_disp_dist", "mean": "estimation"})
beast_data["method"] = "BEAST"

bd_data = bd_data.rename(columns={"max_distance": "max_disp_dist", "calc_diff_rate": "estimation"})
bd_data["method"] = "Brownian Diffusion"

# Split genepop into slope and inv_slope datasets
genepop_slope = genepop_data[["max_distance", "slope", "ci_lower", "ci_upper"]].copy()
genepop_slope = genepop_slope.rename(columns={"max_distance": "max_disp_dist", "slope": "estimation"})
genepop_slope["method"] = "Regression_slope"

genepop_inv = genepop_data[["max_distance", "inv_slope", "inv_ci_lower", "inv_ci_upper"]].copy()
genepop_inv = genepop_inv.rename(columns={"max_distance": "max_disp_dist", "inv_slope": "estimation"})
genepop_inv["method"] = "Regression_invSlope"

# Combine all
combined_df = pd.concat([beast_data, bd_data, genepop_slope, genepop_inv], ignore_index=True)
combined_df = combined_df[["max_disp_dist", "method", "estimation"]]
combined_df = combined_df.sort_values(by=["method", "max_disp_dist"], ignore_index=True)

combined_df.to_csv("combined_data.tsv", index=False, sep="\t")

# Convert max_disp_dist and sample_size to categorical for ANOVA
combined_df["max_disp_dist"] = combined_df["max_disp_dist"].astype("category")

output = []

# 1-Way ANOVA across all methods
output.append("=== 1-Way ANOVA: estimation ~ method (all data) ===")
output.append("# Question: Are the estimation methods significantly different overall across all max dispersal distances?")
method_groups = [g["estimation"].values for _, g in combined_df.groupby("method")]
stat, p = f_oneway(*method_groups)
output.append(f"F = {stat:.4f}, p = {p:.4e}")
output.append("# ✅ Significant difference between methods" if p < 0.05 else "# ❌ No significant difference between methods")
output.append("")

# 1-Way ANOVA: estimation ~ max_disp_dist (per method)
output.append("=== 1-Way ANOVA: estimation ~ max_disp_dist (within each method) ===")
output.append("# Question: Does estimation vary with max dispersal distance within each method?")
for method in combined_df["method"].unique():
    subset = combined_df[combined_df["method"] == method]
    groups = [g["estimation"].values for _, g in subset.groupby("max_disp_dist")]
    if all(len(g) > 1 for g in groups):
        stat, p = f_oneway(*groups)
        output.append(f"Method: {method} --> F = {stat:.4f}, p = {p:.4e}")
        output.append("# ✅ Significant variation within the method" if p < 0.05 else "# ❌ No significant variation within method")
output.append("")

# Linear model: estimation ~ max_disp_dist (per method)
output.append("=== Linear Model: estimation ~ max_disp_dist (within each method) ===")
output.append("# Question: Does a linear trend exist between estimation and max dispersal distance for each method?")
for method in combined_df["method"].unique():
    subset = combined_df[combined_df["method"] == method].copy()
    if subset["max_disp_dist"].dtype.name == "category":
        subset["max_disp_dist"] = subset["max_disp_dist"].astype(float)
    X = sm.add_constant(subset["max_disp_dist"])
    y = subset["estimation"]
    model = sm.OLS(y, X).fit()
    slope = model.params["max_disp_dist"]
    pval = model.pvalues["max_disp_dist"]
    output.append(f"Method: {method} --> slope = {slope:.4f}, p = {pval:.4e}")
    output.append("# ✅ Significant linear trend" if pval < 0.05 else "# ❌ No significant linear trend")
output.append("")

# Save output
with open("stats_for_max_disp_dist.txt", "w") as f:
    f.writelines(line + "\n" for line in output)