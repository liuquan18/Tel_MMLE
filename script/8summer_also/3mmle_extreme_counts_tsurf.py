# %%
import src.MMLE_TEL.story_line as story_line
import numpy as np
import importlib
import matplotlib.pyplot as plt
import src.plots.extrc_tsurf_scatter as extrc_tsurf
import src.MMLE_TEL.index_stats as index_stats
import xarray as xr
import os

#%%
import importlib

importlib.reload(story_line)
importlib.reload(extrc_tsurf)

# %%
import multiprocessing as mp
from multiprocessing import Pool

#%%
# MPI_GE_onepct
# MPI_GE_onepct summer
index_stats.extreme_counts_tsurf(
    model        =  "MPI_GE_onepct", 
    tsurf        =  "ens_fld_year_ocean_mean", 
    standard     =  "first",
    season       =  'MJJA')


# %%
# MPI_GE_onepct winter
index_stats.extreme_counts_tsurf(
    model        = "MPI_GE_onepct", 
    tsurf        = "ens_fld_year_ocean_mean", 
    standard     = "first",
    season       = 'DJFM')

#%%
index_stats.extreme_counts_tsurf(
    model        =  "MPI_GE_onepct", 
    tsurf        =  "ens_fld_year_ocean_mean", 
    standard     =  "temporal_ens",
    season       =  'MJJA')


#%%
index_stats.extreme_counts_tsurf(
    model        = "MPI_GE_onepct", 
    tsurf        = "ens_fld_year_ocean_mean", 
    standard     = "temporal_ens",
    season       = 'DJFM')






#%%
# MMLEA
models = ["MPI_GE_onepct", "MPI_GE", "CanESM2", "CESM1_CAM5", "GFDL_CM3", "MK36"]
tsurfs = ["ens_fld_year_ocean_mean"] #,"ens_fld_year_mean", "NA_tsurf","tropical_arctic_gradient"]
standards = ["first"]

def run_extreme_counts_tsurf(args):
    model, tsurf, standard = args
    index_stats.extreme_counts_tsurf(model, tsurf=tsurf, standard=standard)

with Pool() as p:
    p.map(run_extreme_counts_tsurf, [(model, tsurf,standard) for model in models for tsurf in tsurfs for standard in standards])


