#%%
import src.MMLE_TEL.index_generator as index_generate
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import src.compute.slurm_cluster as scluster
#%%
import importlib
importlib.reload(index_generate)
importlib.reload(scluster)


# %%
# function for generate the index
def index_gen(model,standard = 'first', season = 'MJJA'):
    generator = index_generate.decompose_troposphere(model,standard=standard,season = season)
    generator.save_result()
#%%
index_gen('MPI_GE_onepct',standard = 'first',season = 'MJJA')


# %%
