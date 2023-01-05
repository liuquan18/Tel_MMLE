#!/bin/bash
# activate environment
conda activate thirdPanel

# go to 
cd /work/mh0033/m300883/Tele_season/docs
# upload
swift upload Tel_season build
