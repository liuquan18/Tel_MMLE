# %%
import src.MMLE_TEL.story_line as story_line
import numpy as np
import importlib
import matplotlib.pyplot as plt
import src.plots.extreme_plot as extrc_tsurf
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
    tsurf        =  "ens_fld_year_mean", 
    standard     =  "first",
    season       =  'JJAS')


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
    season       =  'JJAS')


#%%
index_stats.extreme_counts_tsurf(
    model        = "MPI_GE_onepct", 
    tsurf        = "ens_fld_year_ocean_mean", 
    standard     = "temporal_ens",
    season       = 'DJFM')

#%%
# season = 'MAM'
index_stats.extreme_counts_tsurf(
    model        = "MPI_GE_onepct",
    tsurf        = "ens_fld_year_ocean_mean",
    standard     = "first",
    season       = 'MAM')

#%%
# try fixed_pattern = 'all'
index_stats.extreme_counts_tsurf(
    model        = "MPI_GE_onepct",
    tsurf        = "ens_fld_year_ocean_mean",
    fixed_pattern= 'all',
    standard     = "first",
    season       = 'JJA')


#%%
# MMLEA
def run_extreme_counts_tsurf(args):
    model, tsurf, standard, season = args
    index_stats.extreme_counts_tsurf(model, tsurf=tsurf, standard=standard, season=season)

def mmle_extreme_counts_tsurf(season):
    models = ["MPI_GE_onepct", "MPI_GE", "CanESM2","GFDL_CM3", "CESM1_CAM5", "MK36"]
    tsurfs = ["ens_fld_year_mean"] #,"ens_fld_year_mean", "NA_tsurf","tropical_arctic_gradient"]
    standards = ["first"]
    seasons = [season]

    with Pool() as p:
        p.map(run_extreme_counts_tsurf, [(model, tsurf, standard, season) for model in models for tsurf in tsurfs for standard in standards for season in seasons])


# %%
mmle_extreme_counts_tsurf('JJAS')
#%%
mmle_extreme_counts_tsurf('DJFM')


# %%
index_stats.extreme_counts_tsurf(
    model        = "CanESM2",
    tsurf        = "ens_fld_year_mean",
    standard     = "first",
    season       = 'JJAS')
# %%
index_stats.extreme_counts_tsurf(
    model        = "CESM1_CAM5",
    tsurf        = "ens_fld_year_mean",
    standard     = "first",
    season       = 'JJAS')
# %%
index_stats.extreme_counts_tsurf(
    model        = "mk36",
    tsurf        = "ens_fld_year_mean",
    standard     = "first",
    season       = 'JJAS')
# %%
index_stats.extreme_counts_tsurf(
    model        = "MK36",
    tsurf        = "ens_fld_year_mean",
    standard     = "first",
    season       = 'MAM')
# %%
# %%
index_stats.extreme_counts_tsurf(
    model        = "MPI_GE",
    tsurf        = "ens_fld_year_mean",
    standard     = "first",
    season       = 'MAM')
# %%
index_stats.extreme_counts_tsurf(
    model        = "GFDL_CM3",
    tsurf        = "ens_fld_year_mean",
    standard     = "first",
    season       = 'DJFM')
# %%
