#!/bin/bash

module load cdo

FROM=/work/mh0033/m300883/3rdPanel/data/onepct/echam/
TO=/work/mh0033/m300883/3rdPanel/data/influence/t2min/

for ens in {0001..0069}  # not all are avaiable yet.
do
    echo "***************************"
    echo "ens:" ${ens}

	Files=${FROM}onepct${ens}_echam6_echam_1*.grb

	cdo -P  48 -f nc -yearmean -selmonth,1,2,3,4 -selyear,1851/1999 -shifttime,1mo -mergetime -apply,"-selname,var201,var202" [ ${Files} ] ${TO}onepct${ens}_echam6_echam.nc

done