# %%
import src.MMLE_TEL.story_line as story_line
import importlib
import src.plots.extreme_plot as extrc_tsurf
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
    season      = "JJAS",
    )
# %%
index_stats.extreme_counts_profile(
    model       = "MPI_GE_onepct",
    standard    = "first",
    season      = "DJFM",
    )
#%%
index_stats.extreme_counts_profile(
    model       = "MPI_GE_onepct",
    standard    = "temporal_ens",
    season      = "JJAS",
    )

# %%
index_stats.extreme_counts_profile(
    model       = "MPI_GE_onepct",
    standard    = "temporal_ens",
    season      = "DJFM",
    )
# %%
# for season 'MAM'
index_stats.extreme_counts_profile(
    model       = "MPI_GE_onepct",
    standard    = "first",
    season      = "MAM",
    )
# %%
