#%%
import xarray as xr
import numpy as np
import src.MMLE_TEL.story_line as story_line

# %%
mpige_onepct_ind_decade = story_line.story_line("MPI_GE_onepct", "ind", "decade")
# %%
mpige_onepct_ind_decade.plot_all()
# %%