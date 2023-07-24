# %%
# %%
import src.MMLE_TEL.index_generator as index_generate
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import src.compute.slurm_cluster as scluster
import concurrent.futures
import sys
from mpi4py import MPI

# %%
import importlib

importlib.reload(index_generate)
importlib.reload(scluster)


# %%
# function for generate the index
def index_gen(model, fixedPattern = 'decade', plev=50000, standard="first"):
    generator = index_generate.decompose_monthly(
        model, plev=plev, fixedPattern=fixedPattern, standard=standard
    )
    generator.save_result()


# %%
index_gen("MPI_GE_onepct", "decade", 50000, "first")
# %%

eof_result = rolling_eof.rolling_eof(
    g
)