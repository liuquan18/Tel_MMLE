# %%
import src.MMLE_TEL.story_line as story_line
import numpy as np
import importlib
import matplotlib.pyplot as plt
import src.plots.extreme_plot as extrc_tsurf
import src.MMLE_TEL.index_stats as index_stats
import xarray as xr
import os
import sys
#%%
import importlib

importlib.reload(index_stats)
importlib.reload(extrc_tsurf)

# %%
import multiprocessing as mp
from multiprocessing import Pool

#%%
def read_var_months(model):
    odir = f'/work/mh0033/m300883/Tel_MMLE/data/{model}/'
    Jun_fl = odir + 'ts_Jun/'
    Jul_fl = odir + 'ts_Jul/'
    Aug_fl = odir + 'ts_Aug/'

    Jun_f = index_stats.read_var_data(Jun_fl)
    Jul_f = index_stats.read_var_data(Jul_fl)
    Aug_f = index_stats.read_var_data(Aug_fl)

    JJA_f = xr.concat([Jun_f, Jul_f, Aug_f], dim='time')
    return JJA_f

def composite(model):
    var_data = read_var_months(model)
    index_stats.composite_analysis(
        model            = model,
        index_season     = 'JJA',
        tsurf_season     = 'JJA',
        var_data         = var_data,)
    
#%%
num = int(sys.argv[1])
t1 = int(sys.argv[2])
t2 = int(sys.argv[3])



models = ['MPI_GE','CanESM2','CESM1_CAM5','MK36','GFDL_CM3']

# %%
print("===========================================")
print(f"node_num:{num} is doing {models[num-1]}")
composite(models[num-1])

# %%
composite('MPI_GE')
# %%
for model in models:
    composite(model)

# %%
composite('CESM1_CAM5')
# %%
composite('MK36')
# %%
composite('CanESM2')
# %%
composite('GFDL_CM3')
# %%
