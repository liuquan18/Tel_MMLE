# imports 
# %%
import xarray as xr
from src.MMLE_TEL.index_generator import decompose_troposphere
import logging
# %%
model = 'MPI_GE'
# %%
def index_gene(model):
    GEN = decompose_troposphere(
        model = model,
        fixedPattern='decade',
        standard='first',
        season = 'JJA',
        all_years=False,
    )
    
    logging.info(f"Start saving {model}...")
    GEN.save_result()
# %%
index_gene(model)
# %%
