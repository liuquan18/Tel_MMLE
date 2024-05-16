#!/bin/bash
source /etc/profile
module load cdo
module load parallel

export model=$1

# define function zonal_wind
zonal_wind(){
    # input arguments
    infile=$1
    filename=$(basename "$1")

    tmp=/scratch/m/m300883/${model}/zonal_wind_JJA/
    path_out=/work/mh0033/m300883/Tel_MMLE/data/${model}/zonal_wind_JJA/
    mkdir -p $path_out
    mkdir -p $tmp
    # calculate zonal wind
    echo "Calculating zonal wind for $filename"
    cdo -O -fldmean -sellonlatbox,-180,180,20,60 -sellevel,20000 $infile ${tmp}${filename} # zonal wind follow shaw et al, 2023

    # calculate the linear trend of the zonal wind
    echo "Calculating trend for $filename"
    cdo -trend -seasmean -selseason,JJA ${tmp}${filename} ${tmp}a_${filename} ${path_out}trend_${filename}
}

export -f zonal_wind

# a folder to get the name
# for month in Jun, Jul, Aug

Odir=/work/mh0033/m300883/Tel_MMLE/data/${model}/zg/
ls $Odir/*.nc | parallel --jobs 10 zonal_wind {}




