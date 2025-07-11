#!/bin/bash
#SBATCH -c 1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=fast
#SBATCH --mem=8G
#SBATCH -t 03:00:00

source ~/.bashrc
conda activate run_gspace

cd /shared/projects/phylodispsim/results_archive/various_num_of_seqs/500_500_lattice/4seqs
python analyze_reg_results_comp_per_radius.py