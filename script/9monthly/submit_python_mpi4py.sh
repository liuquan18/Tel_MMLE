#!/bin/bash
#SBATCH --job-name=parallel
#SBATCH --time=02:00:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --account=mh0033

mpirun -np 4 python -u hello4mmle_extreme_counts_tsurf.py $1 $2 $3
