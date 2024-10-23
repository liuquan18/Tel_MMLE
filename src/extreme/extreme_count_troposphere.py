#%%
import xarray as xr
import proplot as pplt
import numpy as np


# %%
def extreme_count_rean(model, plevs):
    odir = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/"

    First = []
    Last = []

    for plev in plevs:
        first = xr.open_dataset(f"{odir}first_plev{plev}_extc.nc")
        last = xr.open_dataset(f"{odir}last_plev{plev}_extc.nc")

        # a new dimension called plev
        first['plev'] = plev
        last['plev'] = plev

        First.append(first)
        Last.append(last)


    # concat along the plev dimension
    First = xr.concat(First, dim='plev')
    Last = xr.concat(Last, dim='plev')

    First = First.sortby('plev')
    Last = Last.sortby('plev')

    # divided by ensemble size
    if model == 'CR20_allens':
        first_count = First.pc/(4 * 80)
        last_count = Last.pc/(4 * 80)
    else:
        first_count = First.pc/(4)
        last_count = Last.pc/(4)
    return first_count, last_count

def extreme_count_MPI(model, ens_size=100):
    vertical_eof = 'ind'
    prefix = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/troposphere_{vertical_eof}_decade_first_JJA_"
    first = xr.open_dataset(f"{prefix}first_count.nc")
    last = xr.open_dataset(f"{prefix}last_count.nc")

    # divided by ensemble size
    first_count = first.pc/ens_size
    last_count = last.pc/ens_size
    return first_count, last_count
