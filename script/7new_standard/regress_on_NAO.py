# plot the extreme_count v.s tsurf for all models
# %%
import src.MMLE_TEL.index_stats as index_stats
import numpy as np
import importlib
import xarray as xr
import matplotlib.pyplot as plt
#%%
import importlib
importlib.reload(index_stats)


# %%
def extr_tsurf(model,fixedPattern,ens_size):
    story = index_stats.index_stats(model,vertical_eof = 'ind',fixed_pattern = fixedPattern,standard='temporal_ens',local=True)
    story.to_plot_dir = f"/work/mh0033/m300883/Tel_MMLE/docs/source/plots/MMLE/{model}_NA"
    story.extrc_tsurf(ylim = (0,ens_size/2),ci = 'bootstrap')
    plt.xlim(-1,5)

# %%
extr_tsurf("CanESM2", 'decade', 50)
# %%
extr_tsurf("MPI_GE", 'decade', 100)

# %%
