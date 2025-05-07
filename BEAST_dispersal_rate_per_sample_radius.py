from utils.gspace_utils import *
from utils.beast_utils import *
from utils.misc_utils import *
import random, os, subprocess, sys, shutil

set_tests_dir()

lattice_size_X = int(input("Enter lattice size X: "))
lattice_size_Y = int(input("Enter lattice size Y: "))
radii_to_test = int(input("Enter number of radii to test: "))

radii = [random.randint(1, lattice_size_X//2) for r in range(1, int(radii_to_test)+1)]

for radius in radii:
    # make a dir for each radius
    os.mkdir(f"result_r_{radius}")

    # Generate GSpace settings file per radius specified
    generate_gspace_settings_circular_sample(output_dir=f"./result_r_{radius}",
                                             lattice_size_x=lattice_size_X,
                                             lattice_size_y=lattice_size_Y,
                                             r=radius)

    # change directories for each radius and rename the file
    os.chdir(f"result_r_{radius}")
    os.rename(f"GSpaceSettings_r_{radius}.txt","GSpaceSettings.txt")

    # run GSpace
    run_gspace(gspace_dir="../../../GSpace/build/GSpace")

    # generate the correspoding BEAST XML file
    generate_beast_xml(f"r_{radius}.xml")

    # go back to Tests dir
    os.chdir("..")