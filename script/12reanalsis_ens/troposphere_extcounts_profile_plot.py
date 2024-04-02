#%%
import xarray as xr
import src.plots.extreme_plot as extreme_plot

# %%
odir = "/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/extreme_count/"
# %%
plevs = [92500, 85000, 70000, 50000, 40000, 30000, 20000]
# %%
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
#%%
First = First.sortby('plev')
Last = Last.sortby('plev')
# %%
# divided by ensemble size
first_count = First.pc/(4 * 80)
last_count = Last.pc/(4 * 80)

#%%
extreme_plot.extreme_count_profile(first_count, last_count, xlim = (1,4))

# %%
