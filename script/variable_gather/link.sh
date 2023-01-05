#!/bin/bash

# link the file of onepct to a new directory

FROM=/work/mh1007/MPI-GE/onepct/
TO=/work/mh0033/m300883/3rdPanel/data/onepct/echam/

cd ${TO}

for ens in {0001..0100}
do
    echo "***************************"
    echo "ens:" ${ens}
	Prex=${FROM}onepct${ens}/outdata/echam6/

    for year in {1850..1999}
    do
        Fname=onepct${ens}_echam6_echam_${year}.grb

        Files=${Prex}${Fname}

        # read original path
        Wrongln=$(readlink ${Files})

        # correct the wrong link
        if [[ "$Wrongln" == *"mh1007"* ]]
        then
            Rightln=${Wrongln/mh1007/"mh1007/from_Mistral/mh1007"}

        elif [[ "$Wrongln" == *"mh0033"* ]]
        then
            Rightln=${Wrongln/mh0033/"mh0033/from_Mistral/mh0033"}

        else
            echo "!!!!! not all the above!!!!!"
            echo $Wrongln
            echo $Rightln


        fi

        # link the origal file to the new file
        ln -s ${Rightln} ${Fname}
    done

done

        
    