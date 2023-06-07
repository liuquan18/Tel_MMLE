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
def decompose_random(ens_size,model, fixedPattern = 'decade',plev = 50000):
    index_gen= index_generate.decompose_plev_random_ens(base_model=model,fixedPattern =fixedPattern, ens_size=ens_size,standard='temporal_ens',plev=plev)
    index_gen.save_result()

ens_sizes = [20, 30, 40, 50,60,70,80,90,100]
models = ['MPI_GE_onepct']

#%%
# parallelly implement the decompose_random function with ens_sizes and models
def produce_random_index(plev = 50000):
    with concurrent.futures.ProcessPoolExecutor() as executor:
        futures = []
        for ens_size in ens_sizes:
            for model in models:
                futures.append(executor.submit(decompose_random, ens_size=ens_size,model = model,plev=plev))
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()
            except Exception as exc:
                print(f'generated an exception: {exc}')
# %%
produce_random_index(plev=30000)
# %%
