#%%
import numpy as np
import pandas as pd
import xarray as xr
from tqdm.notebook import tqdm, trange

import src.Teleconnection.spatial_pattern as ssp
import src.Teleconnection.tools as tools
import src.warming_stage.warming_stage as warming_stage
# %%
zg_data = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/zg_processed/allens_zg.nc")
zg_data = zg_data.var156
zg_data = tools.split_ens(zg_data)
# %%
zg_data_500 = zg_data.sel(plev=50000)
# %%
def decompose_single_decade(xarr, nmode=2):
    """decompose a single decade."""
    field = xarr.stack(com=("ens", "time"))

    eof_result = ssp.doeof(field, nmode=nmode, dim="com")

    return eof_result
# %%
eofs = zg_data_500.resample(time="10AS").map(decompose_single_decade, nmode=2)
# %%
for k,g in zg_data_500.resample(time="10AS"):
    print(k)
    print(g)
    break
# %%
geof = decompose_single_decade(g)
# %%
