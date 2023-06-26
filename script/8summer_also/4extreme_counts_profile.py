# %%
import src.MMLE_TEL.story_line as story_line
import numpy as np
import importlib
import matplotlib.pyplot as plt
import src.plots.extrc_tsurf_scatter as extrc_tsurf
import xarray as xr
import os
import src.extreme.extreme_ci as extreme

#%%
import importlib

importlib.reload(story_line)
importlib.reload(extrc_tsurf)


# %%
def extreme_counts_profile(model,standard = 'first',season = 'MJJA'):
    eof_dir = (f"/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/"            
            + "troposphere_ind_decade_"
            + standard
            + "_"
            + season
            + "_eof_result.nc")
    
    eof_result = xr.open_dataset(eof_dir)

    # PCS of the first and last decade
    first_pc = eof_result.pc.isel(time = slice(0,10))
    last_pc = eof_result.pc.isel(time = slice(-10,None))

    print("calculating the extreme event count")
    # extreme counts of first and last 10 decades
    first_count = extreme.extreme_count_xr(first_pc, ci='bootstrap')
    last_count = extreme.extreme_count_xr(last_pc, ci='bootstrap')

    # save the result
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/"
    prefix = f"troposphere_ind_decade_{standard}_{season}_"
    first_count.to_netcdf(odir + prefix + "first_count.nc")
    last_count.to_netcdf(odir + prefix + "last_count.nc")

# %%
extreme_counts_profile("MPI_GE_onepct")
# %%
