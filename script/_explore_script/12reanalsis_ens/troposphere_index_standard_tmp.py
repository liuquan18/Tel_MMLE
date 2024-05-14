#%%
import xarray as xr
import numpy as np
#%%
from src.reanalysis.decompose import standard_period
# %%
# %%
def split_eof (plev):
    odir = "/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/EOF_result/"
    eof_path = odir + f"all_{plev}_eof.nc"
    eof = xr.open_dataset(eof_path)

    first_eof= eof.sel(time = slice('1850', '1889'))
    last_eof = eof.sel(time = slice('1976', '2015'))

    first_eof_std, last_eof_std = standard_period(first_eof, last_eof)

    first_eof_std.to_netcdf(odir + f"first_plev{plev}_eof_std.nc")
    last_eof_std.to_netcdf(odir + f"last_plev{plev}_eof_std.nc")

# %%
plevs = [92500, 85000, 70000, 50000, 40000, 30000, 20000]
for plev in plevs:
    plev = str(plev)
    split_eof(plev)
# %%
