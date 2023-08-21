# %%
import src.MMLE_TEL.index_generator as index_generate
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import sys
from mpi4py import MPI
import src.Teleconnection.rolling_eof as rolling_eof


# %%
import importlib

importlib.reload(index_generate)
importlib.reload(rolling_eof)
#%%
num = int(sys.argv[1])
t1 = int(sys.argv[2])
t2 = int(sys.argv[3])

#%%
models = ['CESM1_CAM5']#['MPI_GE','CanESM2','CESM1_CAM5','MK36','GFDL_CM3']

# %%
def index_gene(model):
    GEN = index_generate.decompose_plev_JJA(
        model = model,
        fixedPattern='decade_mpi',
    )

    GEN.save_result()

# %%
print("===========================================")
print(f"node_num:{num} is doing {models[num-1]}")
index_gene(models[num-1])