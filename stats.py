import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy.stats import f_oneway

# Load data
beast_data = pd.read_csv("/Users/ayoubrayaneaitallaoua/Desktop/IFB_analysis_result_archive/various_max_disp_dist/500_500_lattice/beast_data.tsv", sep="\t")
genepop_data = pd.read_csv("/Users/ayoubrayaneaitallaoua/Desktop/IFB_analysis_result_archive/various_max_disp_dist/500_500_lattice/genepop_results.tsv", sep="\t")
emp_data = pd.read_csv("/Users/ayoubrayaneaitallaoua/Desktop/IFB_analysis_result_archive/various_max_disp_dist/500_500_lattice/empirical_diffusion_all.tsv", sep="\t")

beast_data = beast_data.rename(columns={"disp_dist_max": "max_disp_dist", "mean": "estimation"})
emp_data = emp_data.rename(columns={"max_distance": "max_disp_dist", "calc_diff_rate": "estimation"})
genepop_data = genepop_data.rename(columns={"max_distance": "max_disp_dist", "slope": "estimation"})


beast_data["method"] = "BEAST"
emp_data["method"] = "Empirical"
genepop_data["method"] = "Regression"

genepop_data["simulation"] = genepop_data.groupby("max_disp_dist").cumcount() + 1

# Combine all
combined_df = pd.concat([beast_data, emp_data, genepop_data], ignore_index=True)
combined_df = combined_df[["max_disp_dist", "simulation", "method", "estimation"]]
combined_df = combined_df.sort_values(by=["method", "max_disp_dist", "simulation"], ignore_index=True)

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

# Save output
with open("stats.txt", "w") as f:
    f.writelines(line + "\n" for line in output)