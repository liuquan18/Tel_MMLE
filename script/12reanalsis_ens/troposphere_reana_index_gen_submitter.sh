#!/bin/bash
#SBATCH --job-name=HGT
#SBATCH --time=08:00:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=0
#SBATCH --mail-type=FAIL
#SBATCH --account=mh0033
#SBATCH --output=boot.%j.out


python troposphere_reanalysis_index_gen.py $1 $2