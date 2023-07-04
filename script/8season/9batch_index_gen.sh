#!/bin/bash
#SBATCH -J indexgen          # Specify job name
#SBATCH -p compute          # Use partition shared
#SBATCH -N 1               # Specify number of nodes (1 for serial applications!)
#SBATCH -n 30         # Specify max. number of tasks to be invoked
#SBATCH -t 08:00:00       # Set a limit on the total run time
#SBATCH -A mh0033          # Charge resources on this project account
#SBATCH --mem=0G
#SBATCH -o indexgen.o%j       # File name for standard and error output
#SBATCH --mail-type=FAIL

cd /work/mh0033/m300883/Tel_MMLE/script/8season
source ~/.bashrc
source activate TelSeason

LD_LIBRARY_PATH=/home/m/m300883/libraries
module load openmpi/4.1.2-intel-2021.5.0


srun -n 30 --mpi=pmi2 env MPICC=/sw/spack-levante/openmpi-4.1.2-yfwe6t/bin/mpicc python -m mpi4py hello1index_generator_500hpa.py