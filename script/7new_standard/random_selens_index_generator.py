#%%
#%%
import src.MMLE_TEL.index_generator as index_generate
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
# %%
import importlib
importlib.reload(index_generate)
# %%
import concurrent.futures
# %%
def decompose_random(ens_size,fixedPattern = 'decade'):
    index_gen= index_generate.decompose_plev_random_ens(base_model= 'MPI_GE',fixedPattern =fixedPattern, ens_size=ens_size,standard='temporal_ens')
    index_gen.save_result()

ens_sizes = [20, 30, 40, 50]

with concurrent.futures.ProcessPoolExecutor() as executor:
    futures = [executor.submit(decompose_random, ens_size) for ens_size in ens_sizes]

for future in concurrent.futures.as_completed(futures):
    result = future.result()
# %%
