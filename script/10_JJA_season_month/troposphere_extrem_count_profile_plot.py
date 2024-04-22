
#%%
import xarray as xr
import src.plots.extreme_plot as extreme_plot
# %%
import importlib
importlib.reload(extreme_plot)
# %%
model = 'MPI_GE_onepct'
vertical_eof = 'ind'
prefix = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/troposphere_{vertical_eof}_decade_first_JJA_"
first = xr.open_dataset(f"{prefix}first_count.nc")
last = xr.open_dataset(f"{prefix}last_count.nc")
# %%
# divided by ensemble size
first_count = first.pc/100
last_count = last.pc/100

#%%
extreme_plot.extreme_count_profile(first_count, last_count, xlim = (1,4))

# %%
