#%%
import xarray as xr
import pandas as pd
import numpy as np
# %%
extrc = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/CanESM2/extreme_count/plev_50000_decade_first_JJAS_extre_counts.nc")
tsurf = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/CanESM2/extreme_count/ens_fld_year_mean.nc")
extrc_rands = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_random/extreme_count/plev_50000_decade_JJAS_first_50_extre_counts.nc")
tsurf_rands = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_random/extreme_count/ens_fld_year_mean.nc")


# %%
mode = 'NAO'
extr_type = 'pos'
# %%
import matplotlib.pyplot as plt
import proplot as pplt


#%%

ens_size = 50

def line_single(extrc, tsurf, extrc_rand, tsurf_rand, ens_size, ax, color = 'k'):
    tsurf = tsurf.values
    tsurf = tsurf - tsurf[0]
    extr = extrc.sel(confidence = 'true').pc.squeeze().values / ens_size

    tsurf_rand = tsurf_rand.values
    tsurf_rand = tsurf_rand - tsurf_rand[0]
    extr_rand = extrc_rand.sel(confidence = 'true').pc.squeeze().values / ens_size

    ax.line(x = tsurf,y = extr, marker = 'o', color = color,markersize = 3)
    ax.line(x = tsurf_rand,y = extr_rand, marker = 'o', color = color, alpha = 0.5,markersize = 3)
    ax.set_xlim(-1,5.6)

#%%


extrc = extrc.sel(mode = mode, extr_type = extr_type)
tsurf = tsurf.ts


extrc_rand = extrc_rands.sel(mode = mode, extr_type = extr_type)
tsurf_rands = tsurf_rands.tsurf
#%%

fig = pplt.figure()
ax = fig.subplots()
line_single(extrc, tsurf,extrc_rand,tsurf_rands,  ens_size = 50, ax = ax)
# %%
