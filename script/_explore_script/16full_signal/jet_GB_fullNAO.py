#%%
import xarray as xr
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import os
from src.mechnisms.mechisms import read_greenland_blocking, read_jetStream, Jet_location


# %%
model = 'MPI_GE'
ppc_first_ano = xr.open_dataset(f"/work/mh0033/m300883/Tel_MMLE/data/{model}/full_signal/first_pattern_projected_ano.nc").pc
ppc_last_ano = xr.open_dataset(f"/work/mh0033/m300883/Tel_MMLE/data/{model}/full_signal/last_pattern_projected_ano.nc").pc
#%%

ppc_first_ano = ppc_first_ano.isel(mode = 0)
ppc_last_ano = ppc_last_ano.isel(mode = 0)

#%%
jet_stream = read_jetStream(model)
jet_stream.load()
jet_stream_loc = Jet_location(jet_stream)
#%%
GB = read_greenland_blocking(model)
GB.load()
#%%
def anomaly(x):
    return (x - x.sel(time = slice('1850','1859')).mean())
#%%
jet_stream_loc_ano = anomaly(jet_stream_loc)
GB_ano = anomaly(GB)
#%%
def corresponding_NAO(NAO, var, threshold = 0): 
    """
    var: jet_stream_loc or GB
    """
    std = var.sel(time = slice('1850','1859')).std()
    var_pos = var.where(var > threshold*std)
    NAO_var_pos = NAO.where(var_pos.notnull())

    return NAO_var_pos

#%%
first_NAO_jet_north = corresponding_NAO(ppc_first_ano, jet_stream_loc_ano)
last_NAO_jet_north = corresponding_NAO(ppc_last_ano, jet_stream_loc_ano)

first_NAO_GB_above = corresponding_NAO(ppc_first_ano, GB_ano)
last_NAO_GB_above = corresponding_NAO(ppc_last_ano, GB_ano)

# %%
fig, axes = plt.subplots(3,1, figsize = (5,10))

ppc_first_ano.plot.hist(ax = axes[0], alpha = 0.7,bins=np.arange(-4, 4.1, 0.5),color = 'k', label = 'first10')
ppc_last_ano.plot.hist(ax = axes[0], alpha = 0.7, bins=np.arange(-4, 4.1, 0.5),color = 'r', label = 'last10')

# vline at mean
axes[0].axvline(ppc_first_ano.mean(), color='k', linestyle='dashed', linewidth=1.5)
axes[0].axvline(ppc_last_ano.mean(), color='r', linestyle='dashed', linewidth=1.5)


first_NAO_jet_north.plot.hist(ax = axes[1], alpha = 0.7,bins=np.arange(-4, 4.1, 0.5),color = 'k', label = 'first10')
last_NAO_jet_north.plot.hist(ax = axes[1], alpha = 0.7, bins=np.arange(-4, 4.1, 0.5),color = 'r', label = 'last10')

# vline at mean
axes[1].axvline(first_NAO_jet_north.mean(), color='k', linestyle='dashed', linewidth=1.5)
axes[1].axvline(last_NAO_jet_north.mean(), color='r', linestyle='dashed', linewidth=1.5)


first_NAO_GB_above.plot.hist(ax = axes[2], alpha = 0.7,bins=np.arange(-4, 4.1, 0.5),color = 'k', label = 'first10')
last_NAO_GB_above.plot.hist(ax = axes[2], alpha = 0.7, bins=np.arange(-4, 4.1, 0.5),color = 'r', label = 'last10')

# vline at mean
axes[2].axvline(first_NAO_GB_above.mean(), color='k', linestyle='dashed', linewidth=1.5)
axes[2].axvline(last_NAO_GB_above.mean(), color='r', linestyle='dashed', linewidth=1.5)

for ax in axes.flatten():
    ax.set_xlim(-4.2, 4.2)
    ax.set_ylim(0, 580)

axes[0].set_title("full change in NAO")
axes[1].set_title("NAO with jet stream location more northly than climatology")
axes[2].set_title("NAO with GB index higher than climatology")

axes[0].legend(frameon = False)

# a, b, c
axes[0].text(-5.5, 550, f"{chr(97)}", fontsize=12, fontweight='bold')
axes[1].text(-5.5, 550, f"{chr(98)}", fontsize=12, fontweight='bold')
axes[2].text(-5.5, 550, f"{chr(99)}", fontsize=12, fontweight='bold')

plt.tight_layout()
plt.savefig(f"/work/mh0033/m300883/Tel_MMLE/docs/source/plots/mechism/NAO_full_signal_ano_{model}.pdf")
# %%
