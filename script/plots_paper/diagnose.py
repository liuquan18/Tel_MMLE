#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
# %%
def read_blocking(model):
    path = f'/work/mh0033/m300883/Tel_MMLE/data/{model}/GB_index_JJA'
    blocking = xr.open_mfdataset(f'{path}/trend*.nc', combine = 'nested', concat_dim = 'ens')
    blocking = blocking.squeeze()
    try:
        blocking = blocking.to_dataframe().reset_index()[['zg']]
    except KeyError:
        try:
            blocking = blocking.to_dataframe().reset_index()[['gph']]
        except KeyError:
            blocking = blocking.to_dataframe().reset_index()[['var156']]
    blocking.columns = [model]
    return blocking
# %%
def read_zonal(model):
    path = f'/work/mh0033/m300883/Tel_MMLE/data/{model}/zonal_wind_JJA'
    zonal = xr.open_mfdataset(f'{path}/trend*.nc', combine = 'nested', concat_dim = 'ens')
    zonal = zonal.squeeze()
    try:
        zonal = zonal.to_dataframe().reset_index()[['zg']]
    except KeyError:
        zonal = zonal.to_dataframe().reset_index()[['u']]
    zonal.columns = [model]
    return zonal
# %%
models = ['MPI_GE','CanESM2','GDDL_CM3','MK36'] #'CESM1_CAM5'
# %%
Blockins = [read_blocking(model) for model in models]
# concat all models
Blockins = pd.concat(Blockins, axis = 1)
# %%
for model in models:
    print(model)
    read_blocking(model)
# %%
