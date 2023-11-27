#!/bin/bash
#SBATCH --job-name=GFDL_count
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mail-type=FAIL
#SBATCH --account=mh0033
#SBATCH --output=my_job.%j.out


# mpirun -np 12 python -u hello4mmle_extreme_counts_tsurf.py $1 $2 $3
mpirun -np 12 python -u hello4mmle_extreme_counts_tsurf.py $1 $2 $3
