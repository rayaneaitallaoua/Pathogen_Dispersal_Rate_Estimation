#!/bin/bash -l
#SBATCH -c 1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=fast
#SBATCH --mem=8G
#SBATCH -t 03:00:00

source ~/.bashrc
conda activate run_gspace

cd /shared/ifbstor1/projects/phylodispsim/Pathogen_Dispersal_Rate_Estimation/dispersal_rate_per_maxDist
python analyze_results_per_maxDist.py