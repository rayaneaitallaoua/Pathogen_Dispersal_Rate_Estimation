from utils.gspace_utils import *
from utils.beast_utils import *
from utils.misc_utils import *
import random, os, subprocess, sys, shutil

radii = list(range(1,50))

set_tests_dir()

for radius in radii:
    # make a dir for each radius
    os.mkdir(f"result_r_{radius}")

    # Generate GSpace settings file per radius specified
    generate_gspace_settings_circular_sample(output_dir=f"./result_r_{radius}",
                                             lattice_size_x=100,
                                             lattice_size_y=100,
                                             r=radius,
                                             mutation_rate=1E-6)

    # change directories for each radius and rename the file
    os.chdir(f"result_r_{radius}")
    os.rename(f"GSpaceSettings_r_{radius}.txt","GSpaceSettings.txt")

    # run GSpace
    run_gspace(gspace_dir="../../../GSpace/build/GSpace")

    # go back to Tests dir
    os.chdir("..")