#%%
import xarray as xr

# %%
da = xr.open_dataset("/work/mh0033/m301027/mean_state/ob_climatology/drifter_annualmeans.nc")

# change the longitude values from -180 to 180 to 0 to 360
da = da.sortby('Lon')
lon = da.Lon.values

lon[lon < 0] = lon[lon < 0] + 360

u = xr.DataArray(da.U.T,coords = {'lat':da.Lat.values,'lon':lon},dims = ['lat','lon'])

# %%
