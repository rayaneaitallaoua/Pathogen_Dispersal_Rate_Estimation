import os
import pandas as pd
import matplotlib.pyplot as plt

root = "/Users/ayoubrayaneaitallaoua/Desktop/IFB_analysis_result_archive/various_max_disp_dist/500_500_lattice"

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

if not os.path.exists("plots"):
    os.makedirs("plots")

if not os.path.exists('plots/all_data'):
    os.makedirs('plots/all_data')

if not os.path.exists('plots/inf_20_data'):
    os.makedirs('plots/inf_20_data')

#### ALL DATA ####
### QUANTILES ###
## BEAST - Q ##

# 4seqs #
plt.figure(figsize=(10, 6))
plt.plot(four_seqs_beast['disp_dist_max'], four_seqs_beast["median_diffusion_rate"], marker='o', color='blue', label="4 sequences Median Diffusion Rate")
plt.fill_between(
    four_seqs_beast["disp_dist_max"],
    four_seqs_beast["quantile_2.5"],
    four_seqs_beast["quantile_97.5"],
    color='blue',
    alpha=0.2,
    label="4 sequences 2.5% - 97.5% Quantile Range"
)

# Plot 20seqs #
plt.plot(twenty_seqs_beast['disp_dist_max'], twenty_seqs_beast["median_diffusion_rate"], marker='s', color='green', label="20 sequences Median Diffusion Rate")
plt.fill_between(
    twenty_seqs_beast["disp_dist_max"],
    twenty_seqs_beast["quantile_2.5"],
    twenty_seqs_beast["quantile_97.5"],
    color='green',
    alpha=0.2,
    label="20 sequences 2.5% - 97.5% Quantile Range"
)

# Plot 50seqs
plt.plot(fifty_seqs_beast['disp_dist_max'], fifty_seqs_beast["median_diffusion_rate"], marker='s', color='red', label="50 sequences Median Diffusion Rate")
plt.fill_between(
    fifty_seqs_beast["disp_dist_max"],
    fifty_seqs_beast["quantile_2.5"],
    fifty_seqs_beast["quantile_97.5"],
    color='red',
    alpha=0.2,
    label="50 sequences 2.5% - 97.5% Quantile Range"
)

plt.xlabel("Maximal Dispersal Distance")
plt.ylabel("Median Estimated Diffusion Rate")
plt.title("Estimated Diffusion Rate by Maximal Dispersal Distance for n = (4, 20, 50) [Model: RRW | Method: BEAST]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("plots/all_data/4vs20vs50seqs_beast_quantiles_per_dist_all.png",dpi=300)

## BD Analytical - Q ###

plt.figure(figsize=(10, 6))

# Plot 4seqs
plt.plot(four_seqs_bd['max_distance'], four_seqs_bd["median"], marker='o', color='blue', label="4 sequences Median Diffusion Rate")
plt.fill_between(
    four_seqs_bd["max_distance"],
    four_seqs_bd["quantile_2.5"],
    four_seqs_bd["quantile_97.5"],
    color='blue',
    alpha=0.2,
    label="4 sequences 2.5% - 97.5% Quantile Range"
)

# Plot 20seqs
plt.plot(twenty_seqs_bd['max_distance'], twenty_seqs_bd["median"], marker='s', color='green', label="20 sequences Median Diffusion Rate")
plt.fill_between(
    twenty_seqs_bd["max_distance"],
    twenty_seqs_bd["quantile_2.5"],
    twenty_seqs_bd["quantile_97.5"],
    color='green',
    alpha=0.2,
    label="20 sequences 2.5% - 97.5% Quantile Range"
)

# Plot 50seqs
plt.plot(fifty_seqs_bd['max_distance'], fifty_seqs_bd["median"], marker='s', color='red', label="50 sequences Median Diffusion Rate")
plt.fill_between(
    fifty_seqs_bd["max_distance"],
    fifty_seqs_bd["quantile_2.5"],
    fifty_seqs_bd["quantile_97.5"],
    color='red',
    alpha=0.2,
    label="50 sequences 2.5% - 97.5% Quantile Range"
)

plt.xlabel("Maximal Dispersal Distance")
plt.ylabel("Median Estimated Diffusion Rate")
plt.title("Estimated Diffusion Rate by Maximal Dispersal Distance for n = (4, 20, 50) [Model: Brownian Diffusion | Method: Analytical]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("plots/all_data/4vs20vs50seqs_BD_quantiles_per_dist_all.png",dpi=300)

### IBD REGRESSION - Q ###

plt.figure(figsize=(10, 6))

# Plot 4seqs
plt.plot(four_seqs_ibd_reg['max_distance'], four_seqs_ibd_reg["median"], marker='o', color='blue', label="4 sequences Median Diffusion Rate")
plt.fill_between(
    four_seqs_ibd_reg["max_distance"],
    four_seqs_ibd_reg["quantile_2.5"],
    four_seqs_ibd_reg["quantile_97.5"],
    color='blue',
    alpha=0.2,
    label="4 sequences 2.5% - 97.5% Quantile Range"
)

# Plot 20seqs
plt.plot(twenty_seqs_ibd_reg['max_distance'], twenty_seqs_ibd_reg["median"], marker='s', color='green', label="20 sequences Median Diffusion Rate")
plt.fill_between(
    twenty_seqs_ibd_reg["max_distance"],
    twenty_seqs_ibd_reg["quantile_2.5"],
    twenty_seqs_ibd_reg["quantile_97.5"],
    color='green',
    alpha=0.2,
    label="20 sequences 2.5% - 97.5% Quantile Range"
)

# Plot 50seqs
plt.plot(fifty_seqs_ibd_reg['max_distance'], fifty_seqs_ibd_reg["median"], marker='s', color='red', label="50 sequences Median Diffusion Rate")
plt.fill_between(
    fifty_seqs_ibd_reg["max_distance"],
    fifty_seqs_ibd_reg["quantile_2.5"],
    fifty_seqs_ibd_reg["quantile_97.5"],
    color='red',
    alpha=0.2,
    label="50 sequences 2.5% - 97.5% Quantile Range"
)

plt.xlabel("Maximal Dispersal Distance")
plt.ylabel("Median Estimated Diffusion Rate")
plt.title("Estimated Diffusion Rate by Maximal Dispersal Distance for n = (4, 20, 50) [Model: Isolation By Distance | Method: Regression]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("plots/all_data/4vs20vs50seqs_IBD_reg_quantiles_per_dist_all.png",dpi=300)

### NO QUANTILES ###

## BEAST ##

# 4seqs #
plt.figure(figsize=(10, 6))
plt.plot(four_seqs_beast['disp_dist_max'], four_seqs_beast["median_diffusion_rate"],
         marker='o', color='blue', label="4 sequences Median Diffusion Rate")

# Plot 20seqs #
plt.plot(twenty_seqs_beast['disp_dist_max'], twenty_seqs_beast["median_diffusion_rate"],
         marker='s', color='green', label="20 sequences Median Diffusion Rate")

# Plot 50seqs
plt.plot(fifty_seqs_beast['disp_dist_max'], fifty_seqs_beast["median_diffusion_rate"],
         marker='s', color='red', label="50 sequences Median Diffusion Rate")

plt.xlabel("Maximal Dispersal Distance")
plt.ylabel("Median Estimated Diffusion Rate")
plt.title("Estimated Diffusion Rate by Maximal Dispersal Distance for n = (4, 20, 50) [Model: RRW | Method: BEAST]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("plots/all_data/4vs20vs50seqs_beast_noQuant_per_dist_all.png",dpi=300)

## BD Analytical ###

plt.figure(figsize=(10, 6))

# Plot 4seqs
plt.plot(four_seqs_bd['max_distance'], four_seqs_bd["median"],
         marker='o', color='blue', label="4 sequences Median Diffusion Rate")

# Plot 20seqs
plt.plot(twenty_seqs_bd['max_distance'], twenty_seqs_bd["median"],
         marker='s', color='green', label="20 sequences Median Diffusion Rate")

# Plot 50seqs
plt.plot(fifty_seqs_bd['max_distance'], fifty_seqs_bd["median"],
         marker='s', color='red', label="50 sequences Median Diffusion Rate")

plt.xlabel("Maximal Dispersal Distance")
plt.ylabel("Median Estimated Diffusion Rate")
plt.title("Estimated Diffusion Rate by Maximal Dispersal Distance for n = (4, 20, 50) [Model: Brownian Diffusion | Method: Analytical]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("plots/all_data/4vs20vs50seqs_BD_noQuant_per_dist_all.png",dpi=300)

### IBD REGRESSION ###

plt.figure(figsize=(10, 6))

# Plot 4seqs
plt.plot(four_seqs_ibd_reg['max_distance'], four_seqs_ibd_reg["median"],
         marker='o', color='blue', label="4 sequences Median Diffusion Rate")

# Plot 20seqs
plt.plot(twenty_seqs_ibd_reg['max_distance'], twenty_seqs_ibd_reg["median"],
         marker='s', color='green', label="20 sequences Median Diffusion Rate")

# Plot 50seqs
plt.plot(fifty_seqs_ibd_reg['max_distance'], fifty_seqs_ibd_reg["median"],
         marker='s', color='red', label="50 sequences Median Diffusion Rate")

plt.xlabel("Maximal Dispersal Distance")
plt.ylabel("Median Estimated Diffusion Rate")
plt.title("Estimated Diffusion Rate by Maximal Dispersal Distance for n = (4, 20, 50) [Model: Isolation By Distance | Method: Regression]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("plots/all_data/4vs20vs50seqs_IBD_reg_noQuant_per_dist_all.png",dpi=300)

### BD REGRESSION ###

plt.figure(figsize=(10,6))

# four seqs
plt.plot(four_seqs_bd_reg['disp_dist_max'],four_seqs_bd_reg['slope_median'],
         marker='o', color='blue', label="4 seqs")

# twenty seqs
plt.plot(twenty_seqs_bd_reg['disp_dist_max'],twenty_seqs_bd_reg['slope_median'],
         marker='o', color='green', label="20 seqs")

# fifty seqs
plt.plot(fifty_seqs_bd_reg['disp_dist_max'],fifty_seqs_bd_reg['slope_median'],
         marker='o', color='red', label="50 seqs")

plt.xlabel("Radius")
plt.ylabel("Median Estimated Diffusion Rate")
plt.title("Estimated Diffusion Rate by Maximal Dispersal Distance [Model: Brownian Diffusion | Method: Regression]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("plots/all_data/4vs20vs50seqs_BD_reg_noQuant_per_dist_all.png",dpi=300)

#### inf 20 DATA ####

### filtering and renaming ###

short_four_seqs_beast = four_seqs_beast[four_seqs_beast["disp_dist_max"] < 20]
short_four_seqs_bd = four_seqs_bd[four_seqs_bd["max_distance"] < 20]
short_four_seqs_ibd_reg = four_seqs_ibd_reg[four_seqs_ibd_reg["max_distance"] < 20]
short_four_seqs_bd_reg = four_seqs_bd_reg[four_seqs_bd_reg["disp_dist_max"] < 20]

short_twenty_seqs_beast = twenty_seqs_beast[twenty_seqs_beast["disp_dist_max"] < 20]
short_twenty_seqs_bd = twenty_seqs_bd[twenty_seqs_bd["max_distance"] < 20]
short_twenty_seqs_ibd_reg = twenty_seqs_ibd_reg[twenty_seqs_ibd_reg["max_distance"] < 20]
short_twenty_seqs_bd_reg = twenty_seqs_bd_reg[twenty_seqs_bd_reg["disp_dist_max"] < 20]

short_fifty_seqs_beast = fifty_seqs_beast[fifty_seqs_beast["disp_dist_max"] < 20]
short_fifty_seqs_bd = fifty_seqs_bd[fifty_seqs_bd["max_distance"] < 20]
short_fifty_seqs_ibd_reg = fifty_seqs_ibd_reg[fifty_seqs_ibd_reg["max_distance"] < 20]
short_fifty_seqs_bd_reg = fifty_seqs_bd_reg[fifty_seqs_bd_reg["disp_dist_max"] < 20]

### QUANTILES ###
## BEAST - Q ##

# 4seqs #

plt.figure(figsize=(10, 6))
plt.plot(short_four_seqs_beast['disp_dist_max'], short_four_seqs_beast["median_diffusion_rate"], marker='o', color='blue', label="4 sequences Median Diffusion Rate")
plt.fill_between(
    short_four_seqs_beast["disp_dist_max"],
    short_four_seqs_beast["quantile_2.5"],
    short_four_seqs_beast["quantile_97.5"],
    color='blue',
    alpha=0.2,
    label="4 sequences 2.5% - 97.5% Quantile Range"
)

# Plot 20seqs #
plt.plot(short_twenty_seqs_beast['disp_dist_max'], short_twenty_seqs_beast["median_diffusion_rate"], marker='s', color='green', label="20 sequences Median Diffusion Rate")
plt.fill_between(
    short_twenty_seqs_beast["disp_dist_max"],
    short_twenty_seqs_beast["quantile_2.5"],
    short_twenty_seqs_beast["quantile_97.5"],
    color='green',
    alpha=0.2,
    label="20 sequences 2.5% - 97.5% Quantile Range"
)

# Plot 50seqs
plt.plot(short_fifty_seqs_beast['disp_dist_max'], short_fifty_seqs_beast["median_diffusion_rate"], marker='s', color='red', label="50 sequences Median Diffusion Rate")
plt.fill_between(
    short_fifty_seqs_beast["disp_dist_max"],
    short_fifty_seqs_beast["quantile_2.5"],
    short_fifty_seqs_beast["quantile_97.5"],
    color='red',
    alpha=0.2,
    label="50 sequences 2.5% - 97.5% Quantile Range"
)

plt.xlabel("Maximal Dispersal Distance")
plt.ylabel("Median Estimated Diffusion Rate")
plt.title("Estimated Diffusion Rate by Maximal Dispersal Distance for n = (4, 20, 50) [Model: RRW | Method: BEAST]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("plots/inf_20_data/4vs20vs50seqs_beast_quantiles_per_dist_inf_20.png",dpi=300)

## BD Analytical - Q ###

plt.figure(figsize=(10, 6))

# Plot 4seqs
plt.plot(short_four_seqs_bd['max_distance'], short_four_seqs_bd["median"], marker='o', color='blue', label="4 sequences Median Diffusion Rate")
plt.fill_between(
    short_four_seqs_bd["max_distance"],
    short_four_seqs_bd["quantile_2.5"],
    short_four_seqs_bd["quantile_97.5"],
    color='blue',
    alpha=0.2,
    label="4 sequences 2.5% - 97.5% Quantile Range"
)

# Plot 20seqs
plt.plot(short_twenty_seqs_bd['max_distance'], short_twenty_seqs_bd["median"], marker='s', color='green', label="20 sequences Median Diffusion Rate")
plt.fill_between(
    short_twenty_seqs_bd["max_distance"],
    short_twenty_seqs_bd["quantile_2.5"],
    short_twenty_seqs_bd["quantile_97.5"],
    color='green',
    alpha=0.2,
    label="20 sequences 2.5% - 97.5% Quantile Range"
)

# Plot 50seqs
plt.plot(short_fifty_seqs_bd['max_distance'], short_fifty_seqs_bd["median"], marker='s', color='red', label="50 sequences Median Diffusion Rate")
plt.fill_between(
    short_fifty_seqs_bd["max_distance"],
    short_fifty_seqs_bd["quantile_2.5"],
    short_fifty_seqs_bd["quantile_97.5"],
    color='red',
    alpha=0.2,
    label="50 sequences 2.5% - 97.5% Quantile Range"
)

plt.xlabel("Maximal Dispersal Distance")
plt.ylabel("Median Estimated Diffusion Rate")
plt.title("Estimated Diffusion Rate by Maximal Dispersal Distance for n = (4, 20, 50) [Model: Brownian Diffusion | Method: Analytical]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("plots/inf_20_data/4vs20vs50seqs_BD_quantiles_per_dist_inf_20.png",dpi=300)

### IBD REGRESSION - Q ###

plt.figure(figsize=(10, 6))

# Plot 4seqs
plt.plot(short_four_seqs_ibd_reg['max_distance'], short_four_seqs_ibd_reg["median"], marker='o', color='blue', label="4 sequences Median Diffusion Rate")
plt.fill_between(
    short_four_seqs_ibd_reg["max_distance"],
    short_four_seqs_ibd_reg["quantile_2.5"],
    short_four_seqs_ibd_reg["quantile_97.5"],
    color='blue',
    alpha=0.2,
    label="4 sequences 2.5% - 97.5% Quantile Range"
)

# Plot 20seqs
plt.plot(short_twenty_seqs_ibd_reg['max_distance'], short_twenty_seqs_ibd_reg["median"], marker='s', color='green', label="20 sequences Median Diffusion Rate")
plt.fill_between(
    short_twenty_seqs_ibd_reg["max_distance"],
    short_twenty_seqs_ibd_reg["quantile_2.5"],
    short_twenty_seqs_ibd_reg["quantile_97.5"],
    color='green',
    alpha=0.2,
    label="20 sequences 2.5% - 97.5% Quantile Range"
)

# Plot 50seqs
plt.plot(short_fifty_seqs_ibd_reg['max_distance'], short_fifty_seqs_ibd_reg["median"], marker='s', color='red', label="50 sequences Median Diffusion Rate")
plt.fill_between(
    short_fifty_seqs_ibd_reg["max_distance"],
    short_fifty_seqs_ibd_reg["quantile_2.5"],
    short_fifty_seqs_ibd_reg["quantile_97.5"],
    color='red',
    alpha=0.2,
    label="50 sequences 2.5% - 97.5% Quantile Range"
)

plt.xlabel("Maximal Dispersal Distance")
plt.ylabel("Median Estimated Diffusion Rate")
plt.title("Estimated Diffusion Rate by Maximal Dispersal Distance for n = (4, 20, 50) [Model: Isolation By Distance | Method: Regression]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("plots/inf_20_data/4vs20vs50seqs_IBD_reg_quantiles_per_dist_inf_20.png",dpi=300)

### NO QUANTILES ###

## BEAST ##

# 4seqs #
plt.figure(figsize=(10, 6))
plt.plot(short_four_seqs_beast['disp_dist_max'], short_four_seqs_beast["median_diffusion_rate"],
         marker='o', color='blue', label="4 sequences Median Diffusion Rate")

# Plot 20seqs #
plt.plot(short_twenty_seqs_beast['disp_dist_max'], short_twenty_seqs_beast["median_diffusion_rate"],
         marker='s', color='green', label="20 sequences Median Diffusion Rate")

# Plot 50seqs
plt.plot(short_fifty_seqs_beast['disp_dist_max'], short_fifty_seqs_beast["median_diffusion_rate"],
         marker='s', color='red', label="50 sequences Median Diffusion Rate")

plt.xlabel("Maximal Dispersal Distance")
plt.ylabel("Median Estimated Diffusion Rate")
plt.title("Estimated Diffusion Rate by Maximal Dispersal Distance for n = (4, 20, 50) [Model: RRW | Method: BEAST]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("plots/inf_20_data/4vs20vs50seqs_beast_noQuant_per_dist_inf_20.png",dpi=300)

## BD Analytical ###

plt.figure(figsize=(10, 6))

# Plot 4seqs
plt.plot(short_four_seqs_bd['max_distance'], short_four_seqs_bd["median"],
         marker='o', color='blue', label="4 sequences Median Diffusion Rate")

# Plot 20seqs
plt.plot(short_twenty_seqs_bd['max_distance'], short_twenty_seqs_bd["median"],
         marker='s', color='green', label="20 sequences Median Diffusion Rate")

# Plot 50seqs
plt.plot(short_fifty_seqs_bd['max_distance'], short_fifty_seqs_bd["median"],
         marker='s', color='red', label="50 sequences Median Diffusion Rate")

plt.xlabel("Maximal Dispersal Distance")
plt.ylabel("Median Estimated Diffusion Rate")
plt.title("Estimated Diffusion Rate by Maximal Dispersal Distance for n = (4, 20, 50) [Model: Brownian Diffusion | Method: Analytical]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("plots/inf_20_data/4vs20vs50seqs_BD_noQuant_per_dist_inf_20.png",dpi=300)

### IBD REGRESSION ###

plt.figure(figsize=(10, 6))

# Plot 4seqs
plt.plot(short_four_seqs_ibd_reg['max_distance'], short_four_seqs_ibd_reg["median"],
         marker='o', color='blue', label="4 sequences Median Diffusion Rate")

# Plot 20seqs
plt.plot(short_twenty_seqs_ibd_reg['max_distance'], short_twenty_seqs_ibd_reg["median"],
         marker='s', color='green', label="20 sequences Median Diffusion Rate")

# Plot 50seqs
plt.plot(short_fifty_seqs_ibd_reg['max_distance'], short_fifty_seqs_ibd_reg["median"],
         marker='s', color='red', label="50 sequences Median Diffusion Rate")

plt.xlabel("Maximal Dispersal Distance")
plt.ylabel("Median Estimated Diffusion Rate")
plt.title("Estimated Diffusion Rate by Maximal Dispersal Distance for n = (4, 20, 50) [Model: Isolation By Distance | Method: Regression]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("plots/inf_20_data/4vs20vs50seqs_IBD_reg_noQuant_per_dist_inf_20.png",dpi=300)

### BD REGRESSION ###

plt.figure(figsize=(10,6))

# four seqs
plt.plot(short_four_seqs_bd_reg['disp_dist_max'],short_four_seqs_bd_reg['slope_median'],
         marker='o', color='blue', label="4 seqs")

# twenty seqs
plt.plot(short_twenty_seqs_bd_reg['disp_dist_max'],short_twenty_seqs_bd_reg['slope_median'],
         marker='o', color='green', label="20 seqs")

# fifty seqs
plt.plot(short_fifty_seqs_bd_reg['disp_dist_max'],short_fifty_seqs_bd_reg['slope_median'],
         marker='o', color='red', label="50 seqs")

plt.xlabel("Radius")
plt.ylabel("Median Estimated Diffusion Rate")
plt.title("Estimated Diffusion Rate by Maximal Dispersal Distance [Model: Brownian Diffusion | Method: Regression]")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("plots/inf_20_data/4vs20vs50seqs_BD_reg_noQuant_per_dist_all.png",dpi=300)
