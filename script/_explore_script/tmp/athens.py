#%%
import xarray as xr
import numpy as np
import cartopy.crs as ccrs
import hvplot.xarray
# %%
BE = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE_onepct_30_daily/block_event/onepct_1850-1999_ens_0069.gph500.nc")
# %%
BE.hvplot.contourf(levels=[0.5,1,1.5], geo=True, coastline=True, widget_location='left_top')
# %%
