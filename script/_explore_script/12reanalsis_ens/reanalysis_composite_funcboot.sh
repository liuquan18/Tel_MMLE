#!/bin/bash
#SBATCH --job-name=comp
#SBATCH --time=08:00:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=500G
#SBATCH --mail-type=FAIL
#SBATCH --account=mh0033
#SBATCH --output=comp.%j.out


mpirun -np 1 python -u reananlysis_composite_funcboot.py
# srun -l --cpu_bind=verbose \
#   --distribution=block:cyclic python -u block_event.py $1 $2 $3