# imports 
# %%
import xarray as xr
from src.MMLE_TEL.index_generator import decompose_troposphere
import logging

#%%
import importlib
# %%
# %%
def index_gene(model, vertical_eof = 'dep'):
    GEN = decompose_troposphere(
        model = model,
        fixedPattern='decade',
        standard='first',
        season = 'JJA',
        all_years=False,
        vertical_eof = vertical_eof,
        decades = [1850, 2090]
    )
    
    logging.critical(f"Start saving {model}...")
    GEN.save_result()
# %%
model = 'MPI_GE_onepct'
vertical_eof = 'ind'
index_gene(model, vertical_eof = vertical_eof)
# %%
