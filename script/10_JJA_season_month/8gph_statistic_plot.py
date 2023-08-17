# %%
import xarray as xr
import numpy as np
import os
import sys

import src.MMLE_TEL.index_generator as index_generator
import src.MMLE_TEL.gph_statistic as gph_statistic
import src.Teleconnection.tools as tools



# %%
# %%
models = ["MPI_GE_onepct", "MPI_GE", "CanESM2", "CESM1_CAM5", "MK36", "GFDL_CM3"]

# get the num from keyboard
num = int(sys.argv[1])
# for model in models:
model = models[num - 1]
print(f"calculating the variability of {model} ...")
model_pos, model_neg = gph_statistic.box_variability(model)
to_dir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/box_based/"
# mkdir if not exist
# create directory if it does not exist
if not os.path.exists(to_dir):
    os.makedirs(to_dir)
model_pos.to_netcdf(to_dir + "pos_var.nc")
model_neg.to_netcdf(to_dir + "neg_var.nc")