#!/bin/bash
#SBATCH --job-name=boot
#SBATCH --time=08:00:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=0
#SBATCH --mail-type=FAIL
#SBATCH --account=mh0033
#SBATCH --output=boot.%j.out


mpirun -np 5 python -u reanalysis_composite_boot.py $1 $2 $3
# srun -l --cpu_bind=verbose \
#   --distribution=block:cyclic python -u block_event.py $1 $2 $3