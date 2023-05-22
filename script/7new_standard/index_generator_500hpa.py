#%%
#%%
import xarray as xr
import numpy as np
import pandas as pd

import src.MMLE_TEL.index_generator as index_generate
# %%
MPIGE_decade = index_generate.decompose_mmle("MPI_GE_onepct",plev = 50000,fixedPattern = 'decade')

# %%
MPIGE_decade.save_result()
# %%
MPIGE_all = index_generate.decompose_mmle("MPI_GE_onepct",plev = 50000,fixedPattern = 'all')

# %%
MPIGE_all.save_result()
# %%
