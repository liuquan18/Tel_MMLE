#%%
import xarray as xr
import src.MMLE_TEL.index_generator as index_generate
import numpy as np

# %%
def index_gene(model):
    GEN = index_generate.decompose_plev_JJA(
        model = model,
        fixedPattern='all',
    )

    GEN.save_result()


# %%
index_gene('MPI_GE') # change the src fucntion
# %%
