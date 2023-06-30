#!/bin/bash
#SBATCH -J indexgen          # Specify job name
#SBATCH -p compute          # Use partition shared
#SBATCH -N 1               # Specify number of nodes (1 for serial applications!)
#SBATCH -n 128         # Specify max. number of tasks to be invoked
#SBATCH -t 08:00:00       # Set a limit on the total run time
#SBATCH -A mh0033          # Charge resources on this project account
#SBATCH --mem=100G
#SBATCH -o indexgen.o%j       # File name for standard and error output

cd /work/mh0033/m300883/Tel_MMLE/script/8season
source ~/.bashrc
source activate TelSeason

srun python hello1index_generator_500hpa.py