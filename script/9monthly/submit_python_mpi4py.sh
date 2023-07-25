#!/bin/bash
#SBATCH --job-name=MPI_GE
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --account=mh0033

mpirun -np 12 python -u hello1index_generator_500hpa.py $1 $2 $3
