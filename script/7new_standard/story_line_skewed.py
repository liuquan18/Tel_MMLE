# plot the extreme_count v.s tsurf for all models
# %%
import src.MMLE_TEL.index_stats as index_stats
import numpy as np
import importlib
import matplotlib.pyplot as plt
from scipy.stats import skew

#%%
import importlib
importlib.reload(index_stats)


# %%

import multiprocessing as mp
import src.MMLE_TEL.index_stats as index_stats

def extr_tsurf(model,fixedPattern,ens_size,local = False):
    story = index_stats.index_stats(model,vertical_eof = 'ind',fixed_pattern = fixedPattern,standard='temporal_ens',local = local)
    story.to_plot_dir = f"/work/mh0033/m300883/Tel_MMLE/docs/source/plots/MMLE/{model}_"
    if local:
        story.to_plot_dir = f"/work/mh0033/m300883/Tel_MMLE/docs/source/plots/MMLE/{model}_NA"
    skews = skewness(story.eof_result.pc,mode = 'NAO')
    skews = np.array(skews)
    indices = np.where(skews < 0)
    story.eof_result = story.eof_result.sel(ens = indices[0])
    story.extrc_tsurf(ylim = (0,ens_size/2),ci = 'bootstrap')
    plt.xlim(-1,5)

#%%
def skewness(pc,mode):
    skews = []
    for i in range(0,pc.ens.size):
        skews.append(skew(pc.sel(ens = i,mode = mode)))
    return skews
# %%
extr_tsurf('CanESM2', 'decade', 50)
# %%
extr_tsurf('MPI_GE', 'decade', 50)

# %%