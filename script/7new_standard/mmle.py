# plot the extreme_count v.s tsurf for all models
# %%
import src.MMLE_TEL.index_stats as index_stats
import numpy as np
import importlib
import matplotlib.pyplot as plt

#%%
# same for CanESM2
# Can_decade = index_stats.index_stats("CanESM2",vertical_eof = 'ind',fixed_pattern = 'decade',standard='temporal_ens')
# %%

import multiprocessing as mp
import src.MMLE_TEL.index_stats as index_stats

def extr_tsurf(model,fixedPattern,ens_size):
    story = index_stats.index_stats(model,vertical_eof = 'ind',fixed_pattern = fixedPattern,standard='temporal_ens')
    story.to_plot_dir = f"/work/mh0033/m300883/Tel_MMLE/docs/source/plots/MMLE/{model}_"
    story.extrc_tsurf(ylim = (0,ens_size/2))
    plt.xlim(-1,5)


#%%
if __name__ == '__main__':
    # create a list of arguments for the extr_tsurf function
    args_list = [('CanESM2', 'decade', 50), ('CESM1_CAM5', 'decade', 40), ('MK36', 'decade', 50)]
    
    # create a process for each set of arguments
    processes = [mp.Process(target=extr_tsurf, args=args) for args in args_list]
    
    # start all processes
    for p in processes:
        p.start()
    
    # wait for all processes to finish
    for p in processes:
        p.join()

#%%
extr_tsurf("GFDL_CM3", 'decade', 20)
# %%
