#%%
import xarray as xr
import numpy as np
import pandas as pd

import src.MMLE_TEL.index_generator as index_generate
import src.MMLE_TEL.quick_plot as quick_plot

#%%
# config
v_eof = 'ind' # vertical_eof
fpattern = 'first' # fixed pattern

# %% 
# CANESM2
canesm = index_generate.decompose_fixedPattern("CanESM2",v_eof,fpattern)
canesm.save_result()

# %%
canesm_qp = quick_plot.period_index("CanESM2",v_eof,fpattern,'temp')



# %%