#%%
import xarray as xr
# %%
ds = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/CanESM2/ts/ts_Amon_CanESM2_historical_rcp85_r1i1p1_195001-210012.nc")

# %%
ds.isel(time = 0).ts.plot()
# %%
 