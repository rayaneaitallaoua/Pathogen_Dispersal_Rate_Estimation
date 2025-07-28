import pandas as pd
import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy.stats import f_oneway
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

root = "/Users/ayoubrayaneaitallaoua/Desktop/IFB_analysis_result_archive/various_num_seqs/500_500_lattice"

four_seqs_beast = pd.read_csv(f"{root}/4seqs/beast_diffusion_grouped.tsv", sep="\t")
four_seqs_bd = pd.read_csv(f"{root}/4seqs/empirical_diffusion_grouped.tsv", sep="\t")
four_seqs_bd_reg = pd.read_csv(f"{root}/4seqs/regression_bd_summary.tsv", sep="\t")
four_seqs_genepop_reg_inv = pd.read_csv(f"{root}/4seqs/genepop_invSlope_summary.tsv", sep="\t")

twenty_seqs_beast = pd.read_csv(f"{root}/20seqs/beast_diffusion_grouped.tsv", sep="\t")
twenty_seqs_bd = pd.read_csv(f"{root}/20seqs/empirical_diffusion_grouped.tsv", sep="\t")
twenty_seqs_bd_reg = pd.read_csv(f"{root}/20seqs/regression_bd_summary.tsv", sep="\t")
twenty_seqs_genepop_reg_inv = pd.read_csv(f"{root}/20seqs/genepop_invSlope_summary.tsv", sep="\t")

fifty_seqs_beast = pd.read_csv(f"{root}/50seqs/beast_diffusion_grouped.tsv", sep="\t")
fifty_seqs_bd = pd.read_csv(f"{root}/50seqs/empirical_diffusion_grouped.tsv", sep="\t")
fifty_seqs_bd_reg = pd.read_csv(f"{root}/50seqs/regression_bd_summary.tsv", sep="\t")
fifty_seqs_genepop_reg_inv = pd.read_csv(f"{root}/50seqs/genepop_invSlope_summary.tsv", sep="\t")


# Clean and rename
### 4 SEQS ###
four_seqs_beast = four_seqs_beast.rename(columns={"median": "estimation", "radius": "radius"})
four_seqs_beast["method"] = "BEAST"

four_seqs_bd = four_seqs_bd.rename(columns={"median": "estimation", "radius": "radius"})
four_seqs_bd["method"] = "Brownian Diffusion"

four_seqs_regression_bd = four_seqs_bd_reg[["radius", "slope_median"]].rename(columns={"slope_median": "estimation"})
four_seqs_regression_bd["method"] = "Brownian Diffusion - Regression"

four_seqs_genepop_reg_inv = four_seqs_genepop_reg_inv.rename(columns={"median": "estimation", "radius": "radius"})
four_seqs_genepop_reg_inv["method"] = "IBD regression invSlope (genepop)"

#### 20 SEQS ###

twenty_seqs_beast = twenty_seqs_beast.rename(columns={"median": "estimation", "radius": "radius"})
twenty_seqs_beast["method"] = "BEAST"

twenty_seqs_bd = twenty_seqs_bd.rename(columns={"median": "estimation", "radius": "radius"})
twenty_seqs_bd["method"] = "Brownian Diffusion"

twenty_seqs_regression_bd = twenty_seqs_bd_reg[["radius", "slope_median"]].rename(columns={"slope_median": "estimation"})
twenty_seqs_regression_bd["method"] = "Brownian Diffusion - Regression"

twenty_seqs_genepop_reg_inv = twenty_seqs_genepop_reg_inv.rename(columns={"median": "estimation", "radius": "radius"})
twenty_seqs_genepop_reg_inv["method"] = "IBD regression invSlope (genepop)"

### 50 SEQS ###

fifty_seqs_beast = fifty_seqs_beast.rename(columns={"median": "estimation", "radius": "radius"})
fifty_seqs_beast["method"] = "BEAST"

fifty_seqs_bd = fifty_seqs_bd.rename(columns={"median": "estimation", "radius": "radius"})
fifty_seqs_bd["method"] = "Brownian Diffusion"

fifty_seqs_regression_bd = fifty_seqs_bd_reg[["radius", "slope_median"]].rename(columns={"slope_median": "estimation"})
fifty_seqs_regression_bd["method"] = "Brownian Diffusion - Regression"

fifty_seqs_genepop_reg_inv = fifty_seqs_genepop_reg_inv.rename(columns={"median": "estimation", "radius": "radius"})
fifty_seqs_genepop_reg_inv["method"] = "IBD regression invSlope (genepop)"


output_lines = []

def fit_linear_model(df, method_name):
    df = df.copy()
    df["radius"] = pd.to_numeric(df["radius"])
    X = sm.add_constant(df["radius"])
    y = df["estimation"]
    model = sm.OLS(y, X).fit()
    slope = model.params["radius"]
    pval = model.pvalues["radius"]
    output_lines.append(f"{method_name} — slope: {slope:.5e}, p-value: {pval:.4e}")
    output_lines.append("  ✅ Significant linear trend" if pval < 0.05 else "  ❌ No significant linear trend")

output_lines.append("=== Linear Models: estimation ~ radius ===")

output_lines.append("\n\n*** BEAST ***")
fit_linear_model(four_seqs_beast, "BEAST (4 seqs)")
fit_linear_model(twenty_seqs_beast, "BEAST (20 seqs)")
fit_linear_model(fifty_seqs_beast, "BEAST (50 seqs)")

output_lines.append("\n\n*** Brownian Diffusion ***")
fit_linear_model(four_seqs_bd, "Empirical (4 seqs)")
fit_linear_model(twenty_seqs_bd, "Empirical (20 seqs)")
fit_linear_model(fifty_seqs_bd, "Empirical (50 seqs)")

output_lines.append("\n\n*** Brownian Diffusion (regression) ***")
fit_linear_model(four_seqs_regression_bd, "Brownian Regression (4 seqs)")
fit_linear_model(twenty_seqs_regression_bd, "Brownian Regression (20 seqs)")
fit_linear_model(fifty_seqs_regression_bd, "Brownian Regression (50 seqs)")

output_lines.append("\n\n*** Regression (inverse slope) ***")
fit_linear_model(four_seqs_genepop_reg_inv, "InvRegression (4 seqs)")
fit_linear_model(twenty_seqs_genepop_reg_inv, "InvRegression (20 seqs)")
fit_linear_model(fifty_seqs_genepop_reg_inv, "InvRegression (50 seqs)")


with open("../stats_for_radii.txt", "w") as f:
    f.write("\n".join(output_lines))

#======= FIGURES =======

# Construct slope/p-value summary table from previous computations
slope_data = []

for method_df, label in [
    (four_seqs_beast, "BEAST (4 seqs)"),
    (twenty_seqs_beast, "BEAST (20 seqs)"),
    (fifty_seqs_beast, "BEAST (50 seqs)"),
    (four_seqs_bd, "BD - non phylogenetic (4 seqs)"),
    (twenty_seqs_bd, "BD - non phylogenetic (20 seqs)"),
    (fifty_seqs_bd, "BD - non phylogenetic (50 seqs)"),
    (four_seqs_regression_bd, "BD - regression (4 seqs)"),
    (twenty_seqs_regression_bd, "BD - regression (20 seqs)"),
    (fifty_seqs_regression_bd, "BD - regression (50 seqs)"),
    (four_seqs_genepop_reg_inv, "IBD - regression (4 seqs)"),
    (twenty_seqs_genepop_reg_inv, "IBD - regression (20 seqs)"),
    (fifty_seqs_genepop_reg_inv, "IBD - regression (50 seqs)")
]:
    df = method_df.copy()
    df["radius"] = pd.to_numeric(df["radius"])
    X = sm.add_constant(df["radius"])
    y = df["estimation"]
    model = sm.OLS(y, X).fit()
    slope = model.params["radius"]
    pval = model.pvalues["radius"]
    slope_data.append({"method": label, "slope": slope, "pval": pval})

# Create DataFrame
slope_df = pd.DataFrame(slope_data)
slope_df["-log10(pval)"] = -np.log10(slope_df["pval"])

# Prepare for heatmap
slope_df["type"] = slope_df["method"].str.extract(r"^(.*?)\s\(")
slope_df["seqs"] = slope_df["method"].str.extract(r"\((\d+ seqs)\)")
# Set categorical order for seqs column
slope_df["seqs"] = pd.Categorical(slope_df["seqs"], categories=["4 seqs", "20 seqs", "50 seqs"], ordered=True)

# Pivot slope values
heatmap_data = slope_df.pivot(index="type", columns="seqs", values="slope")

# Create mask for p-values
pval_mask = slope_df.pivot(index="type", columns="seqs", values="pval")
annotations = slope_df.pivot(index="type", columns="seqs", values="slope").applymap(lambda v: f"{v:.3e}")
signif_mask = pval_mask.applymap(lambda p: "**" if p < 0.05 else "")
annotations = annotations + "\n" + signif_mask

# Plot heatmap
plt.figure(figsize=(8, 5))
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
plt.title("Linear slope of estimation ~ radius")
plt.ylabel("Method")
plt.xlabel("Number of sequences")
plt.tight_layout()
plt.savefig("slope_heatmap.png", dpi=300)
plt.close()