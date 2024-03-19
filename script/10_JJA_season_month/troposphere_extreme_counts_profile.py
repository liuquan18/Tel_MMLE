# %%

import xarray as xr
import src.MMLE_TEL.index_stats as index_stats
#%%
import importlib

importlib.reload(index_stats)

# %%
index_stats.extreme_counts_profile(
    model       = "MPI_GE",
    standard    = "first",
    season      = "JJA",
    )
# %%
