#%%
#%%
import src.MMLE_TEL.index_generator as index_generate
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import src.compute.slurm_cluster as scluster
import concurrent.futures

#%%
import importlib
importlib.reload(index_generate)
importlib.reload(scluster)

# %%
# function for generate the index
def index_gen(model,fixedPattern,plev = 50000,season = 'MJJA',standard = 'first'):
    generator = index_generate.decompose_plev(model,plev = plev,fixedPattern = fixedPattern,standard=standard,season = season)
    generator.save_result()
#%%

# for different models
#%%
def index_gen_models(season = 'MAM'):
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [
            # executor.submit(index_gen, 'MPI_GE_onepct', 'decade', plev=50000),
            executor.submit(index_gen, 'MPI_GE', 'decade', plev=50000,season = season,standard = 'first'),
            executor.submit(index_gen, 'CanESM2', 'decade', plev=50000,season = season,standard = 'first'),
            executor.submit(index_gen, 'CESM1_CAM5', 'decade', plev=50000,season = season,standard = 'first'),
            executor.submit(index_gen, 'MK36', 'decade', plev=50000,season = season,standard = 'first'),
            # executor.submit(index_gen, 'GFDL_CM3', 'decade', plev=50000,season = season,standard = 'first'),
        ]
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()
            except Exception as e:
                print(f"Exception: {e}")



# for different seasons and standards, only for MPI_GE_onepct
#%%
def index_gen_seasons(season = 'MAM'):
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = [
            executor.submit(index_gen, 'MPI_GE_onepct', 'decade', plev=50000, season=season, standard='first'),
            executor.submit(index_gen, 'MPI_GE_onepct', 'decade', plev=50000, season=season, standard='first'),
            # executor.submit(index_gen, 'MPI_GE_onepct', 'decade', plev=50000, season='DJFM', standard='first'),
            # executor.submit(index_gen, 'MPI_GE_onepct', 'decade', plev=50000, season='JJAS', standard='temporal_ens'),
            # # executor.submit(index_gen, 'MPI_GE_onepct', 'decade', plev=50000, season='DJFM', standard='temporal_ens')
            # executor.submit(index_gen, 'MPI_GE_onepct', 'decade', plev=50000, season='MAM', standard='temporal_ens')
        ]
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()
            except Exception as e:
                print(f"Exception: {e}")

# define main funtion, run index_gen_seasons() or index_gen_models() here
#%%
def main():
    index_gen_seasons()
    # index_gen_models()
    
# %%

for model in ['MPI_GE','CanESM2','CESM1_CAM5','MK36']:
    print("***************")
    print(model)
    index_gen(model,'decade',plev = 50000,season = 'MAM',standard = 'first')    
# %%