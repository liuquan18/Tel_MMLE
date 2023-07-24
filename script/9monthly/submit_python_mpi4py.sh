#!/bin/bash
#SBATCH --job-name=parallel
#SBATCH --time=02:00:00
#SBATCH --partition=compute
#SBATCH --nodes=2
#SBATCH --ntasks=4
#SBATCH --account=mh0033

mpirun -np 4 python -u hello1index_generator_500hpa.py $1 $2 $3
