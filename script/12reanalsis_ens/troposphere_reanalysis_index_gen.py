# %%
import xarray as xr
import numpy as np
import pandas as pd
import random
import os
import sys
import mpi4py.MPI as MPI

# %%

from src.reanalysis.decompose import EOF_reanalysis
# %%
plev = int(sys.argv[1]) # get the plev from bash file
model = str(sys.argv[2]) # get the model from bash file
# %%
CR20_all_HGT = EOF_reanalysis(model, group_size=40, start_year="1850", end_year="2015", plev = plev, variable = "zg_all_plev")
# %%
CR20_all_HGT.eof.to_netcdf(f"/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/EOF_result/all_{plev}_eof.nc")
CR20_all_HGT.first_eof.to_netcdf(f"/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/EOF_result/first_plev{plev}_eof.nc")
CR20_all_HGT.last_eof.to_netcdf(f"/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/EOF_result/last_plev{plev}_eof.nc")
# %%
