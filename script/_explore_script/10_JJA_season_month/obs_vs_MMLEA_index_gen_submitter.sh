#!/bin/bash
#SBATCH --job-name=obs_vs_MMLEA_index_gen
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=0
#SBATCH --mail-type=FAIL
#SBATCH --account=mh0033
#SBATCH --output=obs_vs_LE.%j.out


mpirun python -u obs_vs_MMLEA_index_gen.py $1
