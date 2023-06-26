# %%
import src.MMLE_TEL.story_line as story_line
import importlib
import src.plots.extrc_tsurf_scatter as extrc_tsurf
import xarray as xr
import src.MMLE_TEL.index_stats as index_stats
#%%
import importlib

importlib.reload(story_line)
importlib.reload(extrc_tsurf)
importlib.reload(index_stats)

# %%
index_stats.extreme_counts_profile(
    model       = "MPI_GE_onepct",
    standard    = "first",
    season      = "MJJA",
    )
# %%
index_stats.extreme_counts_profile(
    model       = "MPI_GE_onepct",
    standard    = "first",
    season      = "DJFM",
    )