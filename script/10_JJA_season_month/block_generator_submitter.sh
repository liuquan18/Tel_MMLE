#!/bin/bash
#SBATCH --job-name=block
#SBATCH --time=08:00:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=0
#SBATCH --mail-type=FAIL
#SBATCH --account=mh0033
#SBATCH --output=block.%j.out


mpirun -np 10 python -u block_generator.py $1 $2 $3