import random

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
Output_Dir=../../TestExample_GSpace/results
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
SampleCoordinateX={sample_x}
SampleCoordinateY={sample_y}
Ind_Per_Node_Sampled={ind_per_node_sampled}
"""
    print("sampled_positions = ")
    print(sampled_positions)
    print("sample_x = ")
    print(sample_x)
    print("sample_y = ")
    print(sample_y)
    print(gspace_settings)

    # Write to file
    with open(f'{output_dir}/GSpaceSettings.txt', "w") as f:
        f.write(gspace_settings)

    print(f"GSpaceSettings.txt generated with random sampling coordinates in {output_dir}!")