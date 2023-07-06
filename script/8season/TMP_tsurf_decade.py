#%%
import xarray as xr
import pandas as pd
import numpy as np
# %%
extrc = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/CanESM2/extreme_count/plev_50000_decade_first_JJAS_extre_counts.nc")
# %%
tsurf = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/CanESM2/ts_processed/ens_fld_year_mean.nc")

# %%
def decade_tsurf(extrc, tsurf, time = 'all'):

    # select time
    if time == 'all':
        pass
    else:
        extrc = extrc.sel(time = slice(time,None))
        tsurf = tsurf.sel(time = slice(time,None))

    time_s = extrc.time.dt.year.values
    time_e = extrc.time[1:].dt.year.values - 1
    time_e = np.append(time_e,2100)

    decade_slices = [slice(str(s), str(e)) for s, e in zip(time_s, time_e)]

    tsurf_dec_mean = [
    tsurf.sel(time=decade_slice).mean(dim="time") for decade_slice in decade_slices
    ]
    tsurf_dec_mean = xr.concat(tsurf_dec_mean, dim=extrc.time)

# %%
