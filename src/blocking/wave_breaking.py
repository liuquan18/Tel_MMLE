# %%
import numpy as np
import xarray as xr

# %%
def wave_breaking_index(Z, lat0_min=30, lat0_max=75, lat_thresh=7.5, lon_thresh=7.5):
    
    WBI = Z.sel(lat=slice( lat0_max,lat0_min)).copy()
    lats = WBI.lat.values
    lons = WBI.lon.values
    
    for i, lat in enumerate(lats):
        lat_s = lat - lat_thresh
        Zlat = Z.sel(lat=lat_s, method='nearest')
        for j, lon in enumerate(lons):
            lon_w = (lon - lon_thresh) % 360
            lon_e = (lon + lon_thresh) % 360
            
            Zw = Zlat.sel(lon=lon_w,method='nearest')
            Ze = Zlat.sel(lon=lon_e,method='nearest')
            wbi = (Zw - Ze) / (2 * lon_thresh)
            WBI.values[:, i, j] = wbi.values
    
    WBI.name = "wave_breaking_index"
    WBI.attrs['long_name'] = "wave breaking index"
    return WBI