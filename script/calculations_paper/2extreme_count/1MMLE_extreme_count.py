# %%
import numpy as np
import importlib
import matplotlib.pyplot as plt
import src.plots.extreme_plot as extrc_tsurf
import src.MMLE_TEL.index_stats as index_stats
import xarray as xr
import os


# %%
import multiprocessing as mp

#%%
# MPI_GE_onepct
# MPI_GE_onepct summer
def extreme_count(model):
    index_stats.extreme_counts(
        model        =  model,
        standard     =  "first",
        season       =  'JJA',
        fixed_pattern=  'all')

# %%
extreme_count('MPI_GE_onepct')
# %%
extreme_count('MPI_GE')
# %%
extreme_count('CanESM2')
# %%
extreme_count('CESM1_CAM5')
# %%
extreme_count('MK36')
# %%
extreme_count('GFDL_CM3')
# %%
