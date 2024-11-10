# %%
import xarray as xr
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from src.mechnisms.mechisms import (
    read_jetStream,
    Jet_location,
    read_greenland_blocking,
    read_NAO_extremes,
)

import matplotlib.pyplot as plt

# %%
jet_stream = read_jetStream("MPI_GE")
jet_stream.load()
jet_loc = Jet_location(jet_stream)

# %%
blocking = read_greenland_blocking("MPI_GE")
blocking.load()

# %%
jet_loc.resample(time = '10Y',closed="left" )
# %%
def decade_corr(x,y):
    x = x.stack(com = ('time','ens'))
    y = y.stack(com = ('time','ens'))
    corr = xr.corr(x, y, dim='com')
    return corr
# %%
jet_GB_corr = jet_loc.resample(time = '10Y',closed="left" ).apply(
    lambda x: decade_corr(x, blocking.sel(time = x.time))
)
# %%
jet_GB_corr = jet_GB_corr.drop_vars(('lon', 'lat', 'plev'))
# %%
fig, ax = plt.subplots()
jet_GB_corr.plot(ax = ax)
ax.set_ylabel("Correlation")
ax.set_title("Jet Stream location and Blocking correlation")
plt.savefig("/work/mh0033/m300883/Tel_MMLE/docs/source/plots/mechism/jet_GB_corr.png")
# %%
