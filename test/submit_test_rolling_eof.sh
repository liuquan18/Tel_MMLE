#!/bin/bash
#SBATCH --job-name=test_rolling
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=10
#SBATCH --account=mh0033

mpirun -np 10 python -u test_rolling_eof.py
