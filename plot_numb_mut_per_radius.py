from utils.gspace_utils import *
from utils.beast_utils import *
from utils.misc_utils import *
import random, os, subprocess, sys, shutil
import pandas as pd
import matplotlib.pyplot as plt

mutation_rate = 1E-5
radii = list(range(1, 50))
results = []

set_tests_dir()

for radius in radii:
    dir_name = f"result_r_{radius}"
    os.mkdir(dir_name)

    # Generate GSpace settings file per radius specified
    generate_gspace_settings_circular_sample(output_dir=f"./{dir_name}",
                                             lattice_size_x=100,
                                             lattice_size_y=100,
                                             r=radius,
                                             mutation_rate=mutation_rate)

    # change directory to run GSpace
    os.chdir(dir_name)
    os.rename(f"GSpaceSettings_r_{radius}.txt", "GSpaceSettings.txt")

    # run GSpace
    run_gspace(gspace_dir="../../../GSpace/build/GSpace")

    # now we are inside result_r_{radius}/, so use local path
    filename = f"sim_seqs_r_{radius}_seq_char_1.txt"
    if os.path.exists(filename):
        df = pd.read_csv(filename, sep="\t")
        avg_mut = df['mut_nb'].mean()
        results.append((radius, avg_mut))
    else:
        print(f"Warning: File not found for radius {radius}")

    os.chdir("..")  # return to Tests dir

# Convert to DataFrame
result_df = pd.DataFrame(results, columns=["radius", "avg_mutations"])
print(result_df)

# Calculate overall average
overall_avg = result_df["avg_mutations"].mean()

# Plot
plt.figure(figsize=(8, 5))
plt.plot(result_df["radius"], result_df["avg_mutations"], marker='o', label='Per-radius average')
plt.axhline(overall_avg, color='red', linestyle='--', label=f'Overall average = {overall_avg:.2f}')
plt.xlabel("Sampling Radius")
plt.ylabel("Average Mutations per Sequence")
plt.title(f"Mutation Load vs Sampling Radius (mutation rate = {mutation_rate})")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig("Number of Mutations vs Sampling Radius.png")