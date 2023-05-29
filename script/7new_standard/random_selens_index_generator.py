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
def standard(generator):
    pc = generator.eof_result['pc']
    std_pc = (pc - pc.mean(dim = ('time','ens')))/pc.std(dim = ('time','ens'))
    generator.std_eof_result['pc'] = std_pc
    return generator


# %%
# random 20 ensemble members
def decompose_random(ens_size,fixedPattern = 'decade'):
    index_gen= index_generate.decompose_plev_random_ens(fixedPattern =fixedPattern, ens_size=ens_size,standard='temporal_ens')
    index_gen.save_result()
# %%
index_20 = decompose_random(20)
# %%
# also for ensembel size of 30, 40, 50
index_30 = decompose_random(30)
index_40 = decompose_random(40)
index_50 = decompose_random(50)
# %%
