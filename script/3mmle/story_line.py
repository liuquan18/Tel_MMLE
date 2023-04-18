#%%
import xarray as xr
import numpy as np
import src.MMLE_TEL.story_line as story_line

import importlib
importlib.reload(story_line)
# %%
mpige_onepct_ind_decade = story_line.story_line("CESM1_CAM5", "ind", "decade")

#%%
mpige_onepct_ind_decade.extreme_count_profile(xlim = (20,120))
#%%
mpige_onepct_ind_decade.plot_all()

#%%
mpige_onepct_ind_decade.create_doc()
# %%
mpige_onepct_dep_decade = story_line.story_line("MPI_GE_onepct", "dep", "decade")
mpige_onepct_dep_decade.extreme_count_profile(xlim = (20,120))
#%%
mpige_onepct_dep_decade.plot_all()
mpige_onepct_dep_decade.create_doc()
# %%
mpige_onepct_dep_decade.stat_overview()
