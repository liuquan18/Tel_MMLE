#%%
import xarray as xr
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import proplot as pplt
import src.plots.extreme_plot as extplt
import src.plots.statistical_overview as stat_overview
import seaborn as sns

#%%
import importlib
importlib.reload(extplt)
importlib.reload(stat_overview)

# %%
def read_extrc(model):
    odir = '/work/mh0033/m300883/Tel_MMLE/data/'+model+'/extreme_count/'
    filename = 'plev_50000_decade_mpi_first_JJA_extre_counts.nc'
    ds = xr.open_dataset(odir+filename).pc
    ds = ds.sel(confidence = 'true')
    return ds
#%%
def read_eof(model):
    odir = f'/work/mh0033/m300883/Tel_MMLE/data/{model}/EOF_result/'
    filename = 'plev_50000_decade_mpi_first_JJA_eof_result.nc'
    ds = xr.open_dataset(odir+filename)
    ds = ds.sel(mode = 'NAO')
    return ds

def split_first_last(eof_result):
    times = eof_result.time
    years = np.unique(times.dt.year)
    first_years = years[:10]
    last_years = years[-10:]

    eof_first = eof_result.isel(decade = 0).sel(time = eof_result['time.year'].isin(first_years))
    eof_last = eof_result.isel(decade = -1).sel(time = eof_result['time.year'].isin(last_years))
    return eof_first,eof_last

# %%
models = ['MPI_GE','CanESM2','CESM1_CAM5','MK36','GFDL_CM3']
extrcs = [read_extrc(model) for model in models]
eof_results = [read_eof(model) for model in models]
eof_firsts,eof_lasts = zip(*[split_first_last(eof_result) for eof_result in eof_results])

#%%
# line plot of extreme counts vs time
extplt.extrc_time_line(extrcs)
plt.savefig(
    '/work/mh0033/m300883/Tel_MMLE/docs/source/plots/monthly/JJA_extreme_counts_line_time_all.png',
)
# %%
# spatial pattern and distribution of index
stat_overview.spatial_index_MMLEA(eof_firsts, eof_lasts)
plt.savefig(
     '/work/mh0033/m300883/Tel_MMLE/docs/source/plots/monthly/JJA_month_spatial_index_stat_MMLEA.png',
)

# %%
# composite analysis
