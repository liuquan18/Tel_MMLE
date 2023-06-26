# %%
import src.MMLE_TEL.story_line as story_line
import numpy as np
import importlib
import matplotlib.pyplot as plt
import src.plots.extrc_tsurf_scatter as extrc_tsurf
import xarray as xr
import os

#%%
import importlib

importlib.reload(story_line)
importlib.reload(extrc_tsurf)


# %%

import multiprocessing as mp


#%%
def extreme_counts_tsurf(model, fixed_pattern="decade", standard = 'temporal_ens',tsurf = 'ens_fld_year_mean',plev = 50000,season = 'MJJA'):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/"
    prefix = f"plev_{plev}_{fixed_pattern}_{standard}_{season}_"
    eof_dir = odir+ "EOF_result/"+ prefix + "eof_result.nc"
    tsurf_dir = odir + "ts_processed/" + tsurf + ".nc"

    # to dir
    extr_counts_dir = odir + "extreme_count/" + prefix + "extre_counts.nc"
    t_surf_mean_dir = odir + "extreme_count/" + tsurf + ".nc"

    # eof
    eof_result = xr.open_dataset(eof_dir)
    # tsurf
    temperature = read_tsurf(tsurf_dir)
    
    # extreme counts
    # check if the extr_counts_dir exists
    try:
        extrc = xr.open_dataset(extr_counts_dir).pc
        _, tsurf_mean = extrc_tsurf.decadal_extrc_tsurf(
            eof_result.pc,ext_counts_xr=extrc, temp = temperature, ci="bootstrap"
        )

    except FileNotFoundError:
        extrc, tsurf_mean = extrc_tsurf.decadal_extrc_tsurf(
            eof_result.pc, temp = temperature, ci="bootstrap"
        )

    # save the result
    try:
        extrc.to_netcdf(extr_counts_dir)
    except PermissionError:
        pass
    tsurf_mean.to_netcdf(t_surf_mean_dir)

    return extrc, tsurf_mean

#%%
# read tsurf
def read_tsurf(tsurf_dir):
    tsurf = xr.open_dataset(tsurf_dir)

    # read the array either tusrf, ts, or tas
    try:
        tsurf_arr = tsurf.tsurf.squeeze()
    except AttributeError:
        try:
            tsurf_arr = tsurf.ts.squeeze()
        except AttributeError:
            tsurf_arr = tsurf.tas.squeeze()
    # change the temp time into datetime64
    try:
        tsurf_arr["time"] = tsurf_arr.indexes["time"].to_datetimeindex()
    except AttributeError:
        pass
    tsurf_arr = tsurf_arr - tsurf_arr[0]
    return tsurf_arr

#%%
from multiprocessing import Pool

models = ["MPI_GE_onepct", "MPI_GE", "CanESM2", "CESM1_CAM5", "GFDL_CM3", "MK36"]
tsurfs = ["ens_fld_year_ocean_mean"] #,"ens_fld_year_mean", "NA_tsurf","tropical_arctic_gradient"]
standards = ["first"]

def run_extreme_counts_tsurf(args):
    model, tsurf, standard = args
    extreme_counts_tsurf(model, tsurf=tsurf, standard=standard)

with Pool() as p:
    p.map(run_extreme_counts_tsurf, [(model, tsurf,standard) for model in models for tsurf in tsurfs for standard in standards])



#%%
# MPI_GE_onepct
extreme_counts_tsurf("MPI_GE_onepct", tsurf="ens_fld_year_ocean_mean", standard="first",season = 'MJJA')


# %%
extreme_counts_tsurf("MPI_GE_onepct", tsurf="ens_fld_year_ocean_mean", standard="first",season = 'DJFM')

# %%
