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
from multiprocessing import Pool

#%%
# MPI_GE_onepct
# MPI_GE_onepct summer
def extreme_count(model):
    index_stats.extreme_counts(
        model        =  model,
        tsurf        =  "ens_fld_year_mean", 
        standard     =  "first",
        season       =  'JJA',
        fixed_pattern=  'decade_mpi')

def tsurf_decade(model):
    index_stats.decadal_tsurf(
        model = model,
        tsurf = "ens_fld_year_mean"
    )

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

#%% tsurf_decade for all models
tsurf_decade('MPI_GE_onepct')
tsurf_decade('MPI_GE')
tsurf_decade('CanESM2')
tsurf_decade('CESM1_CAM5')
tsurf_decade('MK36')
tsurf_decade('GFDL_CM3')

# %%
