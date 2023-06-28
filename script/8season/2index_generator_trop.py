#%%
import src.MMLE_TEL.index_generator as index_generate
import xarray as xr
import numpy as npa
import cartopy.crs as ccrs
import src.compute.slurm_cluster as scluster
#%%
import importlib
importlib.reload(index_generate)
importlib.reload(scluster)


# %%
# function for generate the index
def index_gen(model,standard = 'first', season = 'MJJA',all_years = False):
    generator = index_generate.decompose_troposphere(model,standard=standard,season = season,all_years = all_years)
    generator.save_result()
#%%
index_gen(
    model       = 'MPI_GE_onepct',
    standard    = 'first',
    season      = 'MJJA')


# %%
index_gen(
    model       = 'MPI_GE_onepct',
    standard    = 'first',
    season      = 'DJFM')


#%%
from dask.distributed import Client
client = Client()


futures = []
for season in ['MJJA', 'DJFM']:
    future = client.submit(index_gen, 'MPI_GE_onepct', 'temporal_ens', True, season)
    futures.append(future)

results = client.gather(futures)
# %%
