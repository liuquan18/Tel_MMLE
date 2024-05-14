#%%
import xarray as xr
import numpy as np
import pandas as pd
import cartopy.crs as ccrs

#%%
eof_result = xr.open_dataset(
    "/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct/EOF_result/ind_decade_eof_result.nc"
)
# %%
eof_result.eof
# %%
eof_result.eof.sel(decade="1966", plev=50000).squeeze().plot.contourf(col = 'mode',
    subplot_kws=dict(projection=ccrs.Orthographic(-20,60)),
    levels=np.arange(-1,1.1,0.2),
    transform=ccrs.PlateCarree(),
)
# %%
