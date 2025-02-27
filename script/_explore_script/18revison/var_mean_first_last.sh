#!/bin/bash
#SBATCH --job-name=period
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=200G
#SBATCH --mail-type=FAIL
#SBATCH --account=mh0033
#SBATCH --output=period.%j.out

module load cdo/2.5.0-gcc-11.2.0
var=$1
odir=/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/$var/
todir=/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/${var}_first_last/
files=($odir/*.nc)


# first 1850-1859

# cdo -P 32 -timmean -sellevel,50000 -ensmean -apply,selyear,1850/1859 [ ${files[@]} ] ${todir}${var}_1850-1859.nc

# last 2090-2099

cdo -P 32 -timmean -sellevel,50000 -ensmean -apply,selyear,2090/2099 [ ${files[@]} ] ${todir}${var}_2090-2099.nc