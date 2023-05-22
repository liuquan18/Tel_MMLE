#%%
import xarray as xr
import numpy as np
import pandas as pd

import src.MMLE_TEL.index_generator as index_generate
import src.MMLE_TEL.quick_plot as quick_plot

#%%
import importlib
importlib.reload(index_generate)

#%% 

CESM = index_generate.decompose_mmle("CESM1_CAM5",plev = 50000)
CESM.save_result()
# %%
CANESM = index_generate.decompose_mmle("CanESM2",plev = 50000)
CANESM.save_result()
# %%
MPIGE = index_generate.decompose_mmle("MPI_GE",plev = 50000)
MPIGE.save_result()
# %%
MK36 = index_generate.decompose_mmle("MK36",plev = 50000)
MK36.save_result()
# %%

MPIGE_decade = index_generate.decompose_mmle("MPI_GE",plev = 50000,fixedPattern = 'decade')
# %%
MPIGE_decade.save_result()
# %%
