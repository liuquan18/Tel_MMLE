#!/bin/bash
#SBATCH --job-name=evcount
#SBATCH --time=08:00:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=2
#SBATCH --mem=0
#SBATCH --mail-type=FAIL
#SBATCH --account=mh0033
#SBATCH --output=evcount.%j.out


mpirun -np 5 python -u block_event_counts.py $1 $2 $3
# srun -l --cpu_bind=verbose \
#   --distribution=block:cyclic python -u block_event.py $1 $2 $3