# %%
import src.MMLE_TEL.story_line as story_line
import numpy as np
import importlib
import matplotlib.pyplot as plt
import src.plots.extreme_plot as extrc_tsurf
import src.MMLE_TEL.index_stats as index_stats
import xarray as xr
import os

#%%
import importlib

importlib.reload(index_stats)
importlib.reload(extrc_tsurf)

# %%
import multiprocessing as mp
from multiprocessing import Pool
# %%
index_stats.composite_analysis(
    model            = "MPI_GE_onepct",
    index_season     = 'JJAS',
    tsurf_season     = 'JJAS')
# %%
