#%%
import xarray as xr
import numpy as np
import src.MMLE_TEL.story_line as story_line

import importlib
importlib.reload(story_line)
# %%
mpige_onepct_ind_decade = story_line.story_line("MPI_GE_onepct", "ind", "decade",prefix='plev_50000_all_none_')

#%%
mpige_onepct_ind_decade.extrc_tsurf(ylim=(30, 110))
# %%
mpige_onepct_ind_decade.stat_overview()
# %%
mpige_onepct_ind_all = story_line.story_line("MPI_GE_onepct", "ind", "all")

# %%
