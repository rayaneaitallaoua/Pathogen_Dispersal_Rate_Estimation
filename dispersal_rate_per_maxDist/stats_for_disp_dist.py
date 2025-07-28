import pandas as pd
import statsmodels.api as sm
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

root = "/Users/ayoubrayaneaitallaoua/Desktop/IFB_analysis_result_archive/various_max_disp_dist/500_500_lattice"

#### ALL DATA ####

### Loading ###

four_seqs_beast = pd.read_csv(
    f"{root}/4seqs/beast_data_grouped.tsv", sep="\t")
four_seqs_bd = pd.read_csv(
    f"{root}/4seqs/non_phylo_bd_grouped.tsv", sep="\t")
four_seqs_ibd_reg = pd.read_csv(
    f"{root}/4seqs/genepop_invSlope_summary.tsv", sep="\t")
four_seqs_bd_reg = pd.read_csv(
    f"{root}/4seqs/regression_bd_summary.tsv", sep="\t")

twenty_seqs_beast = pd.read_csv(
    f"{root}/20seqs/beast_data_grouped.tsv", sep="\t")
twenty_seqs_bd = pd.read_csv(
    f"{root}/20seqs/non_phylo_bd_grouped.tsv", sep="\t")
twenty_seqs_ibd_reg = pd.read_csv(
    f"{root}/20seqs/genepop_invSlope_summary.tsv", sep="\t")
twenty_seqs_bd_reg = pd.read_csv(f"{root}/20seqs/regression_bd_summary.tsv", sep="\t")

fifty_seqs_beast = pd.read_csv(
    f"{root}/50seqs/beast_data_grouped.tsv", sep="\t")
fifty_seqs_bd = pd.read_csv(
    "/Users/ayoubrayaneaitallaoua/Desktop/IFB_analysis_result_archive/various_max_disp_dist/500_500_lattice/50seqs/empirical_diffusion_grouped.tsv", sep="\t")
fifty_seqs_ibd_reg = pd.read_csv(
    f"{root}/50seqs/genepop_invSlope_summary.tsv", sep="\t")
fifty_seqs_bd_reg = pd.read_csv(f"{root}/50seqs/regression_bd_summary.tsv", sep="\t")

### Cleaning and renaming ###

### 4 SEQS ###
four_seqs_beast = four_seqs_beast.rename(columns={"median_diffusion_rate": "estimation", "disp_dist_max": "max disp dist"})
four_seqs_beast["method"] = "BEAST"
four_seqs_beast["sample size"] = "4"

four_seqs_bd = four_seqs_bd.rename(columns={"median": "estimation", "max_distance": "max disp dist"})
four_seqs_bd["method"] = "Brownian Diffusion - Analytical"
four_seqs_bd["sample size"] = "4"

four_seqs_regression_bd = four_seqs_bd_reg[["disp_dist_max", "slope_median"]].rename(columns={"slope_median": "estimation", "disp_dist_max": "max disp dist"})
four_seqs_regression_bd["method"] = "Brownian Diffusion - Regression"
four_seqs_regression_bd["sample size"] = "4"

four_seqs_ibd_reg = four_seqs_ibd_reg.rename(columns={"median": "estimation", "max_distance": "max disp dist"})
four_seqs_ibd_reg["method"] = "IBD - regression invSlope"
four_seqs_ibd_reg["sample size"] = "4"

#### 20 SEQS ###

twenty_seqs_beast = twenty_seqs_beast.rename(columns={"median_diffusion_rate": "estimation", "disp_dist_max": "max disp dist"})
twenty_seqs_beast["method"] = "BEAST"
twenty_seqs_beast["sample size"] = "20"

twenty_seqs_bd = twenty_seqs_bd.rename(columns={"median": "estimation", "max_distance": "max disp dist"})
twenty_seqs_bd["method"] = "Brownian Diffusion - Analytical"
twenty_seqs_bd["sample size"] = "20"

twenty_seqs_regression_bd = twenty_seqs_bd_reg[["disp_dist_max", "slope_median"]].rename(columns={"slope_median": "estimation", "disp_dist_max": "max disp dist"})
twenty_seqs_regression_bd["method"] = "Brownian Diffusion - Regression"
twenty_seqs_regression_bd["sample size"] = "20"

twenty_seqs_ibd_reg = twenty_seqs_ibd_reg.rename(columns={"median": "estimation", "max_distance": "max disp dist"})
twenty_seqs_ibd_reg["method"] = "IBD - regression invSlope"
twenty_seqs_ibd_reg["sample size"] = "20"

### 50 SEQS ###

fifty_seqs_beast = fifty_seqs_beast.rename(columns={"median_diffusion_rate": "estimation", "disp_dist_max": "max disp dist"})
fifty_seqs_beast["method"] = "BEAST"
fifty_seqs_beast["sample size"] = "50"

fifty_seqs_bd = fifty_seqs_bd.rename(columns={"median": "estimation", "max_distance": "max disp dist"})
fifty_seqs_bd["method"] = "Brownian Diffusion - Analytical"
fifty_seqs_bd["sample size"] = "50"

fifty_seqs_regression_bd = fifty_seqs_bd_reg[["disp_dist_max", "slope_median"]].rename(columns={"slope_median": "estimation", "disp_dist_max": "max disp dist"})
fifty_seqs_regression_bd["method"] = "Brownian Diffusion - Regression"
fifty_seqs_regression_bd["sample size"] = "50"

fifty_seqs_ibd_reg = fifty_seqs_ibd_reg.rename(columns={"median": "estimation", "max_distance": "max disp dist"})
fifty_seqs_ibd_reg["method"] = "IBD - regression invSlope"
fifty_seqs_ibd_reg["sample size"] = "50"


# Combine all
combined_df = pd.concat([four_seqs_beast, four_seqs_bd, four_seqs_regression_bd, four_seqs_ibd_reg,
                         twenty_seqs_beast, twenty_seqs_bd, twenty_seqs_regression_bd, twenty_seqs_ibd_reg,
                         fifty_seqs_beast, fifty_seqs_bd, fifty_seqs_regression_bd, fifty_seqs_ibd_reg], ignore_index=True)

combined_df = combined_df[["max disp dist", "sample size", "method", "estimation"]]
combined_df = combined_df.sort_values(by=["method", "sample size", "max disp dist"], ignore_index=True)

combined_df.to_csv("combined_data.tsv", index=False, sep="\t")

output_all = []
"""
# 1-Way ANOVA across all methods
output_all.append("=== 1-Way ANOVA: estimation ~ method (all data) ===")
output_all.append("# Question: Are the estimation methods significantly different overall across all max dispersal distances?")
method_groups = [g["estimation"].values for _, g in combined_df.groupby("method")]
stat, p = f_oneway(*method_groups)
output_all.append(f"F = {stat:.4e}, p = {p:.4e}")
output_all.append("# ✅ Significant difference between methods" if p < 0.05 else "# ❌ No significant difference between methods")
output_all.append("")

# 1-Way ANOVA: estimation ~ max_disp_dist (per method)
output_all.append("=== 1-Way ANOVA: estimation ~ max_disp_dist (within each method) ===")
output_all.append("# Question: Does estimation vary with max dispersal distance within each method?")
for method in combined_df["method"].unique():
    subset = combined_df[combined_df["method"] == method]
    groups = [g["estimation"].values for _, g in subset.groupby("max_disp_dist")]
    if all(len(g) > 1 for g in groups):
        stat, p = f_oneway(*groups)
        output_all.append(f"Method: {method} --> F = {stat:.4e}, p = {p:.4e}")
        output_all.append("# ✅ Significant variation within the method" if p < 0.05 else "# ❌ No significant variation within method")
output_all.append("")
"""

# Linear model: estimation ~ max_disp_dist (per method per sample size)
output_all.append("=== Linear Model: estimation ~ max_disp_dist (within each method) ===")

for method in combined_df["method"].unique():
    output_all.append(f"\n#### For method {method} ###")
    for sample_size in combined_df["sample size"].unique():
        output_all.append(f"### For sample size = {sample_size} ###")
        subset = combined_df[combined_df["method"] == method].copy()
        X = sm.add_constant(subset["max disp dist"])
        y = subset["estimation"]
        model = sm.OLS(y, X).fit()
        slope = model.params["max disp dist"]
        pval = model.pvalues["max disp dist"]
        output_all.append(f"Method: {method} --> slope = {slope:.4e}, p = {pval:.4e}")
        output_all.append("# ✅ Significant linear trend" if pval < 0.05 else "# ❌ No significant linear trend")

output_all.append("")

# Save output_all

with open("./stats_for_max_disp_dist_all.txt", "w") as f:
    f.writelines(line + "\n" for line in output_all)


# Construct slope/p-value summary table from combined_df grouped by method and sample size
slope_data = []

for (method, sample_size), group in combined_df.groupby(["method", "sample size"]):
    X = sm.add_constant(group["max disp dist"])
    y = group["estimation"]
    model = sm.OLS(y, X).fit()
    slope = model.params["max disp dist"]
    pval = model.pvalues["max disp dist"]
    slope_data.append({
        "method": method,
        "seqs": f"n = {sample_size}",
        "slope": slope,
        "pval": pval
    })

slop_df = pd.DataFrame(slope_data)
slope_df = pd.DataFrame(slope_data)
slope_df["-log10(pval)"] = -np.log10(slope_df["pval"])

# Enforce order of sample size (columns) as 4, 20, 50
slope_df["seqs"] = pd.Categorical(
    slope_df["seqs"],
    categories=["n = 4", "n = 20", "n = 50"],
    ordered=True
)

# Pivot slope values
heatmap_data = slope_df.pivot(index="method", columns="seqs", values="slope")

# Create annotations from slope and p-value
pval_mask = slope_df.pivot(index="method", columns="seqs", values="pval")
annotations = heatmap_data.applymap(lambda v: f"{v:.3e}")
signif_mask = pval_mask.applymap(lambda p: "**" if p < 0.05 else "")
annotations = annotations + "\n" + signif_mask

# Plot heatmap
plt.figure(figsize=(10, 6))
sns.heatmap(
    heatmap_data,
    annot=annotations,
    fmt="",
    cmap="RdBu_r",
    center=0,
    vmin=-0.05,
    vmax=0.05,
    cbar_kws={"label": "Slope"}
)
plt.title("Linear slope of estimation ~ max dispersal distance")
plt.ylabel("Method")
plt.xlabel("Number of sequences")
plt.tight_layout()
plt.savefig("plots/all_data/slope_heatmap_max_disp_dist.png", dpi=300)
plt.close()

#### inf 20 DATA ####

# --- Heatmap for subset: max disp dist < 20 ---
subset_df = combined_df[combined_df["max disp dist"] < 20].copy()

slope_data = []
for (method, sample_size), group in subset_df.groupby(["method", "sample size"]):
    X = sm.add_constant(group["max disp dist"])
    y = group["estimation"]
    model = sm.OLS(y, X).fit()
    slope = model.params["max disp dist"]
    pval = model.pvalues["max disp dist"]
    slope_data.append({
        "method": method,
        "seqs": f"n = {sample_size}",
        "slope": slope,
        "pval": pval
    })

slope_df = pd.DataFrame(slope_data)
slope_df["-log10(pval)"] = -np.log10(slope_df["pval"])
slope_df["seqs"] = pd.Categorical(
    slope_df["seqs"],
    categories=["n = 4", "n = 20", "n = 50"],
    ordered=True
)

heatmap_data = slope_df.pivot(index="method", columns="seqs", values="slope")
pval_mask = slope_df.pivot(index="method", columns="seqs", values="pval")
annotations = heatmap_data.map(lambda v: f"{v:.3e}")
signif_mask = pval_mask.map(lambda p: "**" if p < 0.05 else "")
annotations = annotations + "\n" + signif_mask

plt.figure(figsize=(10, 6))
sns.heatmap(
    heatmap_data,
    annot=annotations,
    fmt="",
    cmap="RdBu_r",
    center=0,
    vmin=-0.05,
    vmax=0.05,
    cbar_kws={"label": "Slope"}
)
plt.title("Linear slope of estimation ~ max dispersal distance (disp < 20)")
plt.ylabel("Method")
plt.xlabel("Number of sequences")
plt.tight_layout()
plt.savefig("plots/inf_20_data/slope_heatmap_max_disp_dist_lt20.png", dpi=300)
plt.close()