#%%
import src.MMLE_TEL.index_generator as index_generate
import xarray as xr
import numpy as npa
import cartopy.crs as ccrs
import src.compute.slurm_cluster as scluster
import concurrent.futures

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
def index_gen_trop():
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = []
        futures.append(executor.submit(index_gen, 'MPI_GE_onepct', 'first', 'JJAS'))
        futures.append(executor.submit(index_gen, 'MPI_GE_onepct', 'first', 'MAM'))
        futures.append(executor.submit(index_gen, 'MPI_GE_onepct', 'first', 'DJFM'))
        # futures.append(executor.submit(index_gen, 'MPI_GE_onepct', 'temporal_ens', True, 'MAM'))

        for future in concurrent.futures.as_completed(futures):
            result = future.result()

#%%
def main():
    index_gen_trop()
    
# %%
