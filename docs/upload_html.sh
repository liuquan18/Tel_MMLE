#!/bin/bash
# activate environment
conda activate TelSeason
module load py-python-swiftclient

# go to 
cd /work/mh0033/m300883/Tel_MMLE/docs
# upload
swift upload Tel_MMLE build
