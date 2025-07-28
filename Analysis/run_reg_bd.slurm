#!/bin/bash
#SBATCH -c 1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --partition=fast
#SBATCH --mem=8G
#SBATCH -t 03:00:00

source ~/.bashrc
conda activate run_gspace

python reg_bd.py