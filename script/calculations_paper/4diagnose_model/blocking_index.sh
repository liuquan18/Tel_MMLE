#!/bin/bash
source /etc/profile
module load cdo
module load parallel

export model=$1

# define function blocking_index
blocking_index(){
    # input arguments
    infile=$1
    filename=$(basename "$1")

    tmp=/scratch/m/m300883/${model}/GB_index_JJA/
    path_out=/work/mh0033/m300883/Tel_MMLE/data/${model}/GB_index_JJA/
    mkdir -p $path_out
    mkdir -p $tmp
    # calculate blocking index
    echo "Calculating blocking index for $filename"
    cdo -O -fldmean -sellonlatbox,-80,-20,60,80 -sellevel,50000 $infile ${tmp}${filename} # blocking index following here https://psl.noaa.gov/gcos_wgsp/Timeseries/GBI_UL/ 

    # calculate the linear trend of the blocking index
    echo "Calculating trend for $filename"
    cdo -trend -seasmean -selseason,JJA ${tmp}${filename} ${tmp}a_${filename} ${path_out}trend_${filename}
}

export -f blocking_index

# a folder to get the name
# for month in Jun, Jul, Aug

Odir=/work/mh0033/m300883/Tel_MMLE/data/${model}/zg/
ls $Odir/*.nc | parallel --jobs 10 blocking_index {}


# parallel --jobs 10 blocking_index ::: ls $Odir/*.nc ::: {Jun,Jul,Aug}
# ls $Odir/*.nc | parallel --jobs 10 blocking_index {} ::: Jun Jul Aug

