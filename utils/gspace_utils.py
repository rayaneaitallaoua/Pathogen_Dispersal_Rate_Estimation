import random
import os
import subprocess
import math


def generate_gspace_settings(output_dir=".",
                             lattice_size_x=20,
                             lattice_size_y=20,
                             mutation_rate=1E-5,
                             num_sampled_nodes=4,
                             ind_per_node_sampled=5):
    """
    Generates a GSpaceSettings.txt file with random sampling coordinates.

    Parameters:
    - output_dir (str): Name of the output directory where the GSpaceSettings.txt file will be written.
    - lattice_size_x (int): Size of the lattice in X dimension.
    - lattice_size_y (int): Size of the lattice in Y dimension.
    - num_sampled_nodes (int): Number of distinct nodes to sample.
    """

    # Generate unique random coordinates
    sampled_positions = list()
    while len(sampled_positions) < num_sampled_nodes:
        x = random.randint(1, lattice_size_x)
        y = random.randint(1, lattice_size_y)
        sampled_positions.append((x, y))

    # Separate X and Y coordinates
    sample_x = ",".join(str(pos[0]) for pos in sampled_positions)
    sample_y = ",".join(str(pos[1]) for pos in sampled_positions)

    # GSpace settings content
    gspace_settings = f"""%%%%%%%% SIMULATION SETTINGS %%%%%%%%%%%%%%%
Data_filename=simulated_sequences
Run_Number=1

%%%%%%%% OUTPUT FILE FORMAT SETTINGS %%%%%%%
Output_Dir=.
%Coordinate_file=true
Sequence_characteristics_file=true
Fasta=true
Fasta_Single_Line_Seq=true

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

%%%%%%%%%%%%% SAMPLE SETTINGS %%%%%%%%%%%%%%
Sample_Size_X={len(sample_x.split(","))}
Sample_Size_Y={len(sample_y.split(","))}
SampleCoordinateX={sample_x}
SampleCoordinateY={sample_y}
Ind_Per_Node_Sampled={ind_per_node_sampled}
%%%%%% VARIOUS COMPUTATION OPTION S%%%%%%%%%
Diagnostic_Tables = Effective_Dispersal
"""

    # Write to file
    with open(f'{output_dir}/GSpaceSettings.txt', "w") as f:
        f.write(gspace_settings)

    print(f"GSpaceSettings.txt generated with random sampling coordinates in {output_dir}!")


def run_gspace(gspace_dir="../../GSpace/build/GSpace"):
    # 1. Get GSpace executable
    gspace_executable = f"{gspace_dir}"

    # 2. Check if executable exists and run
    if os.path.isfile(gspace_executable):
        subprocess.run([gspace_executable])
    else:
        print(f"Error: {gspace_executable} not found.")


def generate_gspace_settings_square_sample(output_dir=".",
                                             r=3,
                                             lattice_size_x=20,
                                             lattice_size_y=20,
                                             mutation_rate=1E-5,
                                             num_sampled_nodes=4,
                                             ind_per_node_sampled=5):
    # find the sampling area
    lattice_center = (int(lattice_size_x / 2), int(lattice_size_y / 2))

    sampled_positions = list()

    while len(sampled_positions) < num_sampled_nodes:
        x = random.randint(lattice_center[0] - r, lattice_center[0] + r)
        y = random.randint(lattice_center[1] - r, lattice_center[1] + r)
        sampled_positions.append((x, y))

    # Separate X and Y coordinates
    sample_x = ",".join(str(pos[0]) for pos in sampled_positions)
    sample_y = ",".join(str(pos[1]) for pos in sampled_positions)

    gspace_settings = f"""%%%%%%%% SIMULATION SETTINGS %%%%%%%%%%%%%%%
Data_filename=simulated_sequences_r_{r}
Run_Number=1

%%%%%%%% OUTPUT FILE FORMAT SETTINGS %%%%%%%
Output_Dir=.
%Coordinate_file=true
Sequence_characteristics_file=true
Fasta=true
Fasta_Single_Line_Seq=True

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
Diagnostic_Tables = Effective_Dispersal
"""

    # Write to file
    with open(f'{output_dir}/GSpacesettings_r_{r}.txt', "w") as f:
        f.write(gspace_settings)

    print(f"GSpaceSettings_r_{r}.txt generated with square sampling coordinates in {output_dir}!")



def generate_gspace_settings_circular_sample(output_dir=".",
                                             r=3,
                                             lattice_size_x=20,
                                             lattice_size_y=20,
                                             mutation_rate=1E-5,
                                             num_sampled_nodes=4,
                                             ind_per_node_sampled=5):
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

%%%%%%%% OUTPUT FILE FORMAT SETTINGS %%%%%%%
Output_Dir=.
%Coordinate_file=true
Sequence_characteristics_file=true
Fasta=true
Fasta_Single_Line_Seq=True

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
    with open(f'{output_dir}/GSpaceSettings_r_{r}.txt', "w") as f:
        f.write(gspace_settings)

    print(f"GSpaceSettings_r_{r}.txt generated with circular sampling in {output_dir}!")
    return sampled_positions
