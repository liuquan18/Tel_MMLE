#!/bin/bash
#SBATCH --job-name=cesm_remap
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=5
#SBATCH --mem=0
#SBATCH --mail-type=FAIL
#SBATCH --account=mh0033
#SBATCH --output=cesm_rmap.%j.out


mpirun -np 5 python -u index_generator.py $1 $2 $3
# python -u hello0index_generator.py $1 $2 $3

# srun -n 2 --mpi=pmi2 env MPICC=/sw/spack-levante/openmpi-4.1.2-yfwe6t/bin/mpicc python -m mpi4py hello1index_generator_500hpa.py