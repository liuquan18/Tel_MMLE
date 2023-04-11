#%%
import xarray as xr
import numpy as np
import pandas as pd

import src.Teleconnection.season_eof as season_eof
# %%
zg = season_eof.read_data("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/zg/")
# %%
# NAO box -91	30	82	50 | -63	-5	49	31
NAO_box_neg1 = zg.sel(lat=slice(82,50), lon=slice(-91, 30))
NAO_box_pos1 = zg.sel(lat=slice(49,31), lon=slice(-63, 5))
# %%
# EA box
EA_box_neg1 = zg.sel(lat=slice(	71,59), lon=slice(63,99))
EA_box_neg2 = zg.sel(lat=slice(	61,45), lon=slice(-49,-2))
EA_box_pos1 = zg.sel(lat=slice(	34,9), lon=slice(-70,25))
# %%
# define function to get the fldmean, weighted by cos(lat)
def coslat_weighted_mean(box):
    weights = np.cos(np.deg2rad(box.lat))
    box_weighted = box.weighted(weights)
    weighted_mean = box_weighted.mean(("lon", "lat"))
    return weighted_mean
# %%
# NAO = -8.57*10^-5 - 0.015 * NAO_neg1 + 0.017*NAO_pos1
NAO_neg1 = coslat_weighted_mean(NAO_box_neg1)
NAO_pos1 = coslat_weighted_mean(NAO_box_pos1)
NAO = -8.57e-5 - 0.015 * NAO_neg1 + 0.017*NAO_pos1
# %%
# EA= -5.88×10^(-5)-0.0015∙neg1-0.0011∙neg2 +0.0016∙pos1
EA_neg1 = coslat_weighted_mean(EA_box_neg1)
EA_neg2 = coslat_weighted_mean(EA_box_neg2)
EA_pos1 = coslat_weighted_mean(EA_box_pos1)
EA = -5.88e-5 - 0.0015*EA_neg1 - 0.0011*EA_neg2 + 0.0016*EA_pos1
# %%

NAO.name = 'pc'
EA.name = 'pc'

NAO.to_netcdf('/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/BOX_result/NAO.nc')
EA.to_netcdf('/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/BOX_result/EA.nc')