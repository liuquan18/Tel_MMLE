#%%
import xarray as xr
import numpy as np
import pandas as pd

import matplotlib.pyplot as plt
import proplot as pplt
import src.plots.extreme_plot as extplt
import src.plots.statistical_overview as stat_overview
import seaborn as sns

import src.composite.field_composite as field_composite
#%%
import importlib
importlib.reload(extplt)
importlib.reload(stat_overview)
importlib.reload(field_composite)

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
models = ['MPI_GE_onepct','MPI_GE','CanESM2','CESM1_CAM5','MK36','GFDL_CM3']
extrcs = [read_extrc(model) for model in models]
eof_results = [read_eof(model) for model in models]
eof_firsts,eof_lasts = zip(*[split_first_last(eof_result) for eof_result in eof_results])

#%%
# line plot of extreme counts vs time
extplt.extrc_time_line(extrcs,ylim = (25,315))
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
def read_composite_all_model(var_name = 'ts'):
    models = ['MPI_GE','CanESM2','CESM1_CAM5','MK36','GFDL_CM3']
    firsts = {}
    lasts = {}
    
    for model in models:
        first, last = read_composite(model,var_name)
        firsts[model] = first
        lasts[model] = last
    return firsts,lasts

def read_composite(model,var_name):
    odir = '/work/mh0033/m300883/Tel_MMLE/data/'
    first_name = f'plev_50000_decade_mpi_first_JJA_JJA_first_{var_name}_composite.nc'
    last_name = f'plev_50000_decade_mpi_first_JJA_JJA_last_{var_name}_composite.nc'

    first = xr.open_dataset(f"{odir}{model}/composite/{first_name}")
    last = xr.open_dataset(f"{odir}{model}/composite/{last_name}")
    if var_name == 'ts':
        try:
            first = first.tsurf
            last = last.tsurf
        except AttributeError:
            first = first.ts
            last = last.ts
    elif var_name == 'pr':
        try:
            first = first.pr
            last = last.pr
        except AttributeError:
            first = first.precip
            last = last.precip
    return first,last

# %%
firsts,lasts = read_composite_all_model(var_name='ts')
# %%
field_composite.composite_plot_MMLEA(firsts,lasts,extr_type='pos',levels =np.arange(-1.5,1.6,0.3))
plt.savefig(
    '/work/mh0033/m300883/Tel_MMLE/docs/source/plots/monthly/JJA_month_composite_MMLEA_pos_ts.png',
)
# %%
field_composite.composite_plot_MMLEA(firsts,lasts,extr_type='neg',levels =np.arange(-1.5,1.6,0.3))
plt.savefig(
    '/work/mh0033/m300883/Tel_MMLE/docs/source/plots/monthly/JJA_month_composite_MMLEA_neg_ts.png',
)


# %%
firsts,lasts = read_composite_all_model(var_name='pr')
models = ['MPI_GE','CanESM2','CESM1_CAM5','MK36','GFDL_CM3']

## change unit from mm/s to mm/day
for model in models:
    firsts[model] = firsts[model]*86400
    lasts[model] = lasts[model]*86400
# %%
field_composite.composite_plot_MMLEA(firsts,lasts,extr_type='pos',levels =np.arange(-1.0,1.1,0.2))
plt.savefig(
    '/work/mh0033/m300883/Tel_MMLE/docs/source/plots/monthly/JJA_month_composite_MMLEA_pos_pr.png',
)
# %%
field_composite.composite_plot_MMLEA(firsts,lasts,extr_type='neg',levels =np.arange(-1.0,1.1,0.2))
plt.savefig(
    '/work/mh0033/m300883/Tel_MMLE/docs/source/plots/monthly/JJA_month_composite_MMLEA_neg_pr.png',
)

# %%
