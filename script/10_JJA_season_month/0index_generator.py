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
# %%
MK36 = index_generate.decompose_plev_JJA(
    model = 'MK36',
    fixedPattern='decade_mpi',
)


# %%
