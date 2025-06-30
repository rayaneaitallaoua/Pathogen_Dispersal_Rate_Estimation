import glob
import math
import os
import random
import subprocess

from Bio import SeqIO

mcmc_chain_length = 10000000
g_mut_rate = 1E-6
#radii = [10,20,30,40,50]
radii = [10]
rep_per_radius = 1
nodes_sampled = 20
ind_per_node = 5
latt_size_X = 100
latt_size_Y = 100

def generate_gspace_settings_circular_sample(output_dir=".",
                                             r=3,
                                             lattice_size_x=20,
                                             lattice_size_y=20,
                                             mutation_rate=1E-6,
                                             num_sampled_nodes=4,
                                             ind_per_node_sampled=5,
                                             seed=3):
    # find lattice center
    lattice_center = (lattice_size_x // 2, lattice_size_y // 2)
    sampled_positions = set()

    # find the neighbour node presenting the least distance
    def closest_lattice_node(x, y):

        # candidates are the 4 neighbouring nodes
        candidates = [(math.floor(x), math.floor(y)),
                      (math.floor(x), math.ceil(y)),
                      (math.ceil(x), math.floor(y)),
                      (math.ceil(x), math.ceil(y))]

        # define a lambda  function that returns the distance between our point pt and the candidates,
        # iterates over candidates
        return min(candidates, key=lambda pt: math.dist((x, y), pt))

    # generate random phi and r' values
    while len(sampled_positions) < num_sampled_nodes:
        phi = random.uniform(0, 2 * math.pi)
        r_prime = random.uniform(0, r)

        # convert Polar coordinates to Cartesian
        x = lattice_center[0] + r_prime * math.cos(phi)
        y = lattice_center[1] + r_prime * math.sin(phi)

        node = closest_lattice_node(x, y)

        # Ensure node is inside the lattice bounds
        if 0 <= node[0] < lattice_size_x and 0 <= node[1] < lattice_size_y:
            sampled_positions.add(node)

    sample_x = ",".join(str(pos[0]) for pos in sampled_positions)
    sample_y = ",".join(str(pos[1]) for pos in sampled_positions)

    gspace_settings = f"""%%%%%%%% SIMULATION SETTINGS %%%%%%%%%%%%%%%
Data_filename=sim_seqs_r_{r}
Run_Number=1
Random_seeds={r + seed}

%%%%%%%% OUTPUT FILE FORMAT SETTINGS %%%%%%%
Output_Dir=.
Sequence_characteristics_file=true
Fasta=true
Fasta_Single_Line_Seq=True
Genepop=true

%%%%%%%% MARKERS SETTINGS %%%%%%%%%%%%%%%%%%
Ploidy=Haploid
Chromosome_number=1
Sequence_Size=1000
Mutation_Model=HKY
Mutation_Rate={mutation_rate}

%%%%%%%% RECOMBINATION SETTINGS %%%%%%%%%%%%
Recombination_Rate=0

%%%%%%%% DEMOGRAPHIC SETTINGS %%%%%%%%%%%%%%
%% LATTICE
Lattice_Size_X={lattice_size_x}
Lattice_Size_Y={lattice_size_y}
Ind_Per_Pop=30

%% DISPERSAL
Dispersal_Distribution=uniform
Disp_Dist_Max=1,1
Total_Emigration_Rate=0.05

%%%%%%%% SAMPLE SETTINGS %%%%%%%%%%%%%%%%%%%
Sample_Size_X={len(sample_x.split(","))}
Sample_Size_Y={len(sample_y.split(","))}

%% Lattice center: {lattice_center} %%
SampleCoordinateX={sample_x}
SampleCoordinateY={sample_y}
Ind_Per_Node_Sampled={ind_per_node_sampled}
%%%%%% VARIOUS COMPUTATION OPTION S%%%%%%%%%
%Diagnostic_Tables = Effective_Dispersal%
"""

    # Write to file
    with open(f'{output_dir}/GSpaceSettings.txt', "w") as f:
        f.write(gspace_settings)

    print(f"GSpaceSettings_r_{r}.txt generated with circular sampling in {output_dir}!")
    return sampled_positions

def submit_gspace_slurm_job(gspace_dir="../../GSpace/build/GSpace", job_name="gspace_job", time="00:10:00"):
    slurm_script = f"""#!/bin/bash
#SBATCH -c 1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=fast
#SBATCH --mem=4G
#SBATCH -t 03:00:00
#SBATCH -J {job_name}
#SBATCH -o {job_name}.%j.out
#SBATCH -e {job_name}.%j.err

srun {gspace_dir}
"""

    with open("run_gspace.slurm", "w") as f:
        f.write(slurm_script)

    subprocess.run(["sbatch", "run_gspace.slurm"])

def run_analysis(radii_list, repetitions=10):

    for radius in radii_list:
        for sim in range(repetitions):
            print(f"Radius {radius}, Simulation {sim+1}")
            dirname = f"./r{radius}_sim{sim+1}"
            os.makedirs(dirname, exist_ok=True)
            os.chdir(dirname)

            generate_gspace_settings_circular_sample(lattice_size_x=latt_size_X,
                                                     lattice_size_y=latt_size_Y,
                                                     mutation_rate=g_mut_rate,
                                                     r=radius,
                                                     seed=sim+1,
                                                     num_sampled_nodes=nodes_sampled,
                                                     ind_per_node_sampled=ind_per_node)
            submit_gspace_slurm_job(gspace_dir="../../../GSpace/build/GSpace", job_name=f"gspace_r_{radius}_sim_{sim}", time="00:10:00")

            os.chdir("..")

run_analysis(radii_list=radii,repetitions=rep_per_radius)