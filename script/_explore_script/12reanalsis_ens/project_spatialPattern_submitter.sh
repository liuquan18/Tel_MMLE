#!/bin/bash
#SBATCH --job-name=project
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=0
#SBATCH --mail-type=FAIL
#SBATCH --account=mh0033
#SBATCH --output=project.%j.out


mpirun python -u project_spatialPattern.py $1
