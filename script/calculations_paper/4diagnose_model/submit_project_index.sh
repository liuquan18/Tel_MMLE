#!/bin/bash
#SBATCH --job-name=project
#SBATCH --time=08:00:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --mem=0
#SBATCH --mail-type=FAIL
#SBATCH --account=mh0033
#SBATCH --output=project.%j.out
source ~./bashrc
module load python3/unstable
conda activate /home/m/m300883/miniconda3/envs/TelSeason
mpirun python -u project_index_single_model.py $1 $2 # node-index, model-name