import pandas as pd
import matplotlib.pyplot as plt

root = "."

four_seqs_beast = pd.read_csv(
    f"{root}/4seqs/beast_diffusion_grouped.tsv", sep="\t")
four_seqs_emp = pd.read_csv(
    f"{root}/4seqs/non_phylo_bd_grouped.tsv", sep="\t")
four_seqs_reg = pd.read_csv(
    f"{root}/4seqs/genepop_invSlope_summary.tsv", sep="\t")
four_bd_reg = pd.read_csv(
    f"{root}/4seqs/regression_bd_summary.tsv", sep="\t")


twenty_seqs_beast = pd.read_csv(
    f"{root}/20seqs/beast_diffusion_grouped.tsv", sep="\t")
twenty_seqs_emp = pd.read_csv(
    f"{root}/20seqs/non_phylo_bd_grouped.tsv", sep="\t")
twenty_seqs_reg = pd.read_csv(
    f"{root}/20seqs/genepop_invSlope_summary.tsv", sep="\t")
twenty_bd_reg = pd.read_csv(f"{root}/20seqs/regression_bd_summary.tsv", sep="\t")

fifty_seqs_beast = pd.read_csv(
    f"{root}/50seqs/beast_diffusion_grouped.tsv", sep="\t")
fifty_seqs_emp = pd.read_csv(
    f"{root}/50seqs/non_phylo_bd_grouped.tsv", sep="\t")
fifty_seqs_reg = pd.read_csv(
    f"{root}/50seqs/genepop_invSlope_summary.tsv", sep="\t")
fifty_bd_reg = pd.read_csv(f"{root}/50seqs/regression_bd_summary.tsv", sep="\t")



plt.figure(figsize=(10, 6))

### BEAST - QUANTILES ###

# Plot 4seqs
plt.plot(four_seqs_beast['radius'], four_seqs_beast["median"], marker='o', color='blue', label="4 sequences Median Diffusion Rate")
plt.fill_between(
    four_seqs_beast["radius"],
    four_seqs_beast["quantile_2.5"],
    four_seqs_beast["quantile_97.5"],
    color='blue',
    alpha=0.2,
    label="4 sequences 2.5% - 97.5% Quantile Range"
)

# Plot 20seqs
plt.plot(twenty_seqs_beast['radius'], twenty_seqs_beast["median"], marker='s', color='green', label="20 sequences Median Diffusion Rate")
plt.fill_between(
    twenty_seqs_beast["radius"],
    twenty_seqs_beast["quantile_2.5"],
    twenty_seqs_beast["quantile_97.5"],
    color='green',
    alpha=0.2,
    label="20 sequences 2.5% - 97.5% Quantile Range"
)

# Plot 50seqs
plt.plot(fifty_seqs_beast['radius'], fifty_seqs_beast["median"], marker='s', color='red', label="50 sequences Median Diffusion Rate")
plt.fill_between(
    fifty_seqs_beast["radius"],
    fifty_seqs_beast["quantile_2.5"],
    fifty_seqs_beast["quantile_97.5"],
    color='red',
    alpha=0.2,
    label="50 sequences 2.5% - 97.5% Quantile Range"
)

plt.xlabel("Sampling radii")
plt.ylabel("Median Estimated Diffusion Rate")
plt.title("Estimated Diffusion Rate by Sampling Radius for n = (4, 20, 50) [Model: RRW | Method: BEAST]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("4vs20vs50seqs_beast_quantiles_radii.png",dpi=300)

### BD - QUANTILES ###

plt.figure(figsize=(10, 6))

# Plot 4seqs
plt.plot(four_seqs_emp['radius'], four_seqs_emp["median"], marker='o', color='blue', label="4 sequences Median Diffusion Rate")
plt.fill_between(
    four_seqs_emp["radius"],
    four_seqs_emp["quantile_2.5"],
    four_seqs_emp["quantile_97.5"],
    color='blue',
    alpha=0.2,
    label="4 sequences 2.5% - 97.5% Quantile Range"
)

# Plot 20seqs
plt.plot(twenty_seqs_emp['radius'], twenty_seqs_emp["median"], marker='s', color='green', label="20 sequences Median Diffusion Rate")
plt.fill_between(
    twenty_seqs_emp["radius"],
    twenty_seqs_emp["quantile_2.5"],
    twenty_seqs_emp["quantile_97.5"],
    color='green',
    alpha=0.2,
    label="20 sequences 2.5% - 97.5% Quantile Range"
)

# Plot 50seqs
plt.plot(fifty_seqs_emp['radius'], fifty_seqs_emp["median"], marker='s', color='red', label="50 sequences Median Diffusion Rate")
plt.fill_between(
    fifty_seqs_emp["radius"],
    fifty_seqs_emp["quantile_2.5"],
    fifty_seqs_emp["quantile_97.5"],
    color='red',
    alpha=0.2,
    label="50 sequences 2.5% - 97.5% Quantile Range"
)

plt.xlabel("Sampling radii")
plt.ylabel("Median Estimated Diffusion Rate")
plt.title("Estimated Diffusion Rate by Sampling Radius for  [Model: Brownian Diffusion | Method: Analytical]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("4vs20vs50seqs_BD_quantiles_radii.png",dpi=300)

###### BD REGRESSION - QUANTILES #######

plt.figure(figsize=(10,6))

# four seqs
plt.plot(four_bd_reg['radius'],four_bd_reg['slope_median'],
         marker='o', color='blue', label="4 seqs")
plt.fill_between(
    four_bd_reg["radius"],
    four_bd_reg["slope_q2.5"],
    four_bd_reg["slope_q97.5"],
    color='blue',
    alpha=0.2,
    label="Brownian Diffusion 2.5% - 97.5% Quantile Range"
)

# twenty seqs
plt.plot(twenty_bd_reg['radius'],twenty_bd_reg['slope_median'],
         marker='o', color='green', label="20 seqs")
plt.fill_between(
    twenty_bd_reg["radius"],
    twenty_bd_reg["slope_q2.5"],
    twenty_bd_reg["slope_q97.5"],
    color='green',
    alpha=0.2,
    label="Brownian Diffusion 2.5% - 97.5% Quantile Range"
)

# fifty seqs
plt.plot(fifty_bd_reg['radius'],fifty_bd_reg['slope_median'],
         marker='o', color='red', label="50 seqs")
plt.fill_between(
    fifty_bd_reg["radius"],
    fifty_bd_reg["slope_q2.5"],
    fifty_bd_reg["slope_q97.5"],
    color='red',
    alpha=0.2,
    label="Brownian Diffusion 2.5% - 97.5% Quantile Range"
)

plt.xlabel("Radius")
plt.ylabel("Median Estimated Diffusion Rate")
plt.title("Estimated Diffusion Rate by radius [Model: Brownian Diffusion | Method: Regression]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("4vs20vs50seqs_BD_reg_Quantiles_radii.png",dpi=300)

### IBD REGRESSION - QUANTILES ###

plt.figure(figsize=(10, 6))

# Plot 4seqs
plt.plot(four_seqs_reg['radius'], four_seqs_reg["median"], marker='o', color='blue', label="4 sequences Median Diffusion Rate")
plt.fill_between(
    four_seqs_reg["radius"],
    four_seqs_reg["quantile_2.5"],
    four_seqs_reg["quantile_97.5"],
    color='blue',
    alpha=0.2,
    label="4 sequences 2.5% - 97.5% Quantile Range"
)

# Plot 20seqs
plt.plot(twenty_seqs_reg['radius'], twenty_seqs_reg["median"], marker='s', color='green', label="20 sequences Median Diffusion Rate")
plt.fill_between(
    twenty_seqs_reg["radius"],
    twenty_seqs_reg["quantile_2.5"],
    twenty_seqs_reg["quantile_97.5"],
    color='green',
    alpha=0.2,
    label="20 sequences 2.5% - 97.5% Quantile Range"
)

# Plot 50seqs
plt.plot(fifty_seqs_reg['radius'], fifty_seqs_reg["median"], marker='s', color='red', label="50 sequences Median Diffusion Rate")
plt.fill_between(
    fifty_seqs_reg["radius"],
    fifty_seqs_reg["quantile_2.5"],
    fifty_seqs_reg["quantile_97.5"],
    color='red',
    alpha=0.2,
    label="50 sequences 2.5% - 97.5% Quantile Range"
)

plt.xlabel("Sampling radii")
plt.ylabel("Median Estimated Diffusion Rate")
plt.title("Estimated Diffusion Rate by Sampling Radius for n = (4, 20, 50) [Model: Isolation By Distance | Method: Regression]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("4vs20vs50seqs_IBD_reg_quantiles_radii.png",dpi=300)

######## NO QUANTILES ########

### IBD REGRESSION ###

plt.figure(figsize=(10, 6))

# Plot 4seqs
plt.plot(four_seqs_reg['radius'], four_seqs_reg["median"], marker='o', color='blue', label="4 sequences Median Diffusion Rate")


# Plot 20seqs
plt.plot(twenty_seqs_reg['radius'], twenty_seqs_reg["median"], marker='s', color='green', label="20 sequences Median Diffusion Rate")


# Plot 50seqs
plt.plot(fifty_seqs_reg['radius'], fifty_seqs_reg["median"], marker='s', color='red', label="50 sequences Median Diffusion Rate")

plt.xlabel("Sampling radii")
plt.ylabel("Median Estimated Diffusion Rate")
plt.title("Estimated Diffusion Rate by Sampling Radius for n = (4, 20, 50) [Model: Isolation By Distance | Method: Regression]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("4vs20vs50seqs_IBD_reg_noQuant_radii.png",dpi=300)

###### BD REGRESSION #######

plt.figure(figsize=(10,6))

# four seqs
plt.plot(four_bd_reg['radius'],four_bd_reg['slope_median'],
         marker='o', color='blue', label="4 seqs")

# twenty seqs
plt.plot(twenty_bd_reg['radius'],twenty_bd_reg['slope_median'],
         marker='o', color='green', label="20 seqs")

# fifty seqs
plt.plot(fifty_bd_reg['radius'],fifty_bd_reg['slope_median'],
         marker='o', color='red', label="50 seqs")

plt.xlabel("Radius")
plt.ylabel("Median Estimated Diffusion Rate")
plt.title("Estimated Diffusion Rate by radius [Model: Brownian Diffusion | Method: Regression]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("4vs20vs50seqs_BD_reg_noQuant_radii.png",dpi=300)

