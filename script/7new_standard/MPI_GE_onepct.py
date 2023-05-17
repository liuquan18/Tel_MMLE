#%%
import xarray as xr
import numpy as np
import pandas as pd

import src.MMLE_TEL.index_generator as index_generate
import src.MMLE_TEL.quick_plot as quick_plot

#%%
import importlib
importlib.reload(index_generate)
# %%

index_gen = index_generate.decompose_fixedPattern("MPI_GE_onepct",'ind','decade',standard='none')
index_gen.save_result()
# %%

index_gen = index_generate.decompose_fixedPattern("MPI_GE_onepct",'ind','all',standard='none')
index_gen.save_result()