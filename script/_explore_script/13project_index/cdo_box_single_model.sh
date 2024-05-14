#!/bin/bash

# southern box (90W-40E, 20N-50N)
# northern box (90W-40E, 50N-80N)

# This is a little bit different from McKenna, et.al, 2021. which is (90W-60E, 20N-55N), (90W-60E, 55N-90N)
# This is because the domain where the EOF is applied to calculate the indirect response is (90W-40E, 20N-80N)
source /etc/profile

module load cdo
module load parallel


model=$1
export Odir=/work/mh0033/m300883/Tel_MMLE/data/$model/zg_JJA
export Southdir=/scratch/m/m300883/box/$model/southern_box
export Northdir=/scratch/m/m300883/box/$model/northern_box
export BOXdiffdir=/scratch/m/m300883/box/$model/box_diff
export TO=/work/mh0033/m300883/Tel_MMLE/data/$model/box_diff

mkdir -p $Southdir
mkdir -p $Northdir
mkdir -p $BOXdiffdir
mkdir -p $TO


box_NAO() {
    # 90W-40E, 20N-50N
    infile=$1
    southname=$(basename "$infile" .nc)_south.nc
    northname=$(basename "$infile" .nc)_north.nc

    # calculate the different between the two boxes
    cdo -O -fldmean -sellonlatbox,-90,40,20,50 -sellevel,50000 $infile $Southdir/$southname
    cdo -O -fldmean -sellonlatbox,-90,40,50,80 -sellevel,50000 $infile $Northdir/$northname
    cdo -O -sub $Southdir/$southname $Northdir/$northname $BOXdiffdir/$(basename "$infile" .nc)_box_diff.nc

    # calculate anomaly based on mean and std of first year
    cdo -sub $TO/$(basename "$infile" .nc)_box_diff.nc -timmean -seltimestep,1/10 $BOXdiffdir/$(basename "$infile" .nc)_box_diff.nc $BOXdiffdir/$(basename "$infile" .nc)_sub.nc
    cdo -div $BOXdiffdir/$(basename "$infile" .nc)_sub.nc -timstd -seltimestep,1/10 $BOXdiffdir/$(basename "$infile" .nc)_box_diff.nc $TO/$(basename "$infile" .nc)_anom.nc
    
}

export -f box_NAO

ls $Odir/*.nc | parallel --jobs 5 box_NAO