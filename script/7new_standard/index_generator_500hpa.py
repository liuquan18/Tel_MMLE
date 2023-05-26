#%%
#%%
import xarray as xr
import numpy as np
import pandas as pd

import src.MMLE_TEL.index_generator as index_generate
# %%
MPIGE_decade = index_generate.decompose_plev("MPI_GE_onepct",plev = 50000,fixedPattern = 'decade',standard='temporal_ens')

# %%
MPIGE_decade.save_result()
# %%
MPIGE_all = index_generate.decompose_plev("MPI_GE_onepct",plev = 50000,fixedPattern = 'all',standard='temporal_ens')

# %%
MPIGE_all.save_result()
# %%
CESM_all = index_generate.decompose_plev("CESM1_CAM5",plev = 50000,fixedPattern = 'all',standard='temporal_ens')
# %%
CESM_all.save_result()

# %%
Can_all = index_generate.decompose_plev("CanESM2",plev = 50000,fixedPattern = 'all',standard='temporal_ens')
# %%
Can_all.save_result()
# %%
