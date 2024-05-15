#!/bin/bash
source /etc/profile


# define function blocking_index
blocking_index(){
    # input arguments
    model=$1
    month=$2
    # define path
    path_in=/work/mh0033/m300883/Tel_MMLE/data/${model}/zg_${month}
    path_out=${path}/${model}_${var}_${season}_${year}_${month}_${day}_blocking_index.nc
    # calculate blocking index
    cdo -s -b 32 -O -L setctomiss,0 -gtc,0 -setmissval,0 -setgrid,${path_in} -run,blocking_index,${path_in} ${path_out}
}

odir='/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/zg_Aug