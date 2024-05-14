#%%
import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
# %%
CR20_eof_first = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/EOF_result/first_40_eof_std.nc")
CR20_eof_last = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/CR20_allens/EOF_result/last_40_eof_std.nc")
# %%
fig, axes = plt.subplots(1, 2, figsize=(10, 4))
CR20_eof_first.pc.sel(mode = 'NAO').plot.line(x = 'time',color = 'grey',add_legend = False,ax = axes[0])
CR20_eof_last.pc.sel(mode = 'NAO').plot.line(x = 'time',color = 'grey',add_legend = False,ax = axes[1])
# %%
MPI_GE_eof = xr.open_dataset("/work/mh0033/m300883/Tel_MMLE/data/MPI_GE/EOF_result/plev_50000_decade_mpi_first_JJA_eof_result.nc")
# %%
MPI_GE_pc = MPI_GE_eof.pc.sel(mode = 'NAO')
# %%
MPI_GE_pc = MPI_GE_pc.sel(time = slice('1850','2015'))
# %%
