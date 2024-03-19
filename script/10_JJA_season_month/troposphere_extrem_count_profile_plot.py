
#%%
import xarray as xr
import src.plots.extreme_plot as extreme_plot
# %%
import importlib
importlib.reload(extreme_plot)
# %%
model = 'MPI_GE'
prefix = f"/work/mh0033/m300883/Tel_MMLE/data/{model}/extreme_count/troposphere_ind_decade_first_JJA_"
first = xr.open_dataset(f"{prefix}first_count.nc")
last = xr.open_dataset(f"{prefix}last_count.nc")
# %%
# divided by ensemble size
first_count = first.pc/100
last_count = last.pc/100

#%%
extreme_plot.extreme_count_profile(first_count, last_count, xlim = (0,2))
# %%


# 500 hpa
plev = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/extreme_count/plev_50000_decade_mpi_first_JJA_extre_counts.nc")
# %%
plev_index = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/EOF_result/plev_50000_decade_mpi_first_JJA_eof_result.nc")
# %%
first_index = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/EOF_result/troposphere_ind_decade_first_JJA_eof_result.nc")
# %%
